// p1tool - Phase One IIQ developer used by PlotExtractor (IIQ tabs).
//
// Adapted from APPN-UQ PhaseOne_Pipeline (IIQ_Proc/IIQ2TIFF.cpp, IIQ_Proc/IIQHRclip.cpp),
// Copyright (c) 2025 APPN-UQ, MIT License - see THIRD_PARTY_NOTICES.md.
// Uses the Phase One ImageSDK, which is licensed separately by Phase One A/S and is NOT
// distributed with PlotExtractor: the SDK folder (ImageSdkCpp.dll + SensorProfiles) is
// supplied by the user at run time.
//
// Differences from HR_clip.exe: IIQ files are only READ and can live anywhere (read-only
// shares included); every output path, the settings file and the job list are given
// explicitly, so nothing is ever written next to the IIQ files. One job per line, so any
// number of crops per frame.
//
// Usage:  p1tool <sensor_profiles_dir> <jobs.tsv>
// jobs.tsv, one job per line, TAB separated:
//   mode  iiq_path  out_path  config_json  x  y  w  h  scale  bits  jpeg_quality  compress
//   mode   : tiff | jpg | raw | info
//   x y w h: crop in full-size output pixels (w = 0 -> whole frame)
//   scale  : output scale (1 = full size)
//   bits   : 8 or 16 (16 for tiff/raw only)
//   compress: 0 = none, 1 = LZW (tiff)
// "raw" writes the bare RGB pixels row by row (no header) for PlotExtractor to read.
// Progress lines on stdout start with "@@" (the SDK prints its own timing lines):
//   @@OK <i> <w> <h> <bits>     @@INFO <i> <w> <h> <rawW> <rawH>     @@ERR <i> <message>

#include <P1Image.hpp>
#include <windows.h>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static std::wstring widen(const std::string& s)
{
    if (s.empty()) return std::wstring();
    int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), nullptr, 0);
    std::wstring w(n, 0);
    MultiByteToWideChar(CP_UTF8, 0, s.data(), (int)s.size(), &w[0], n);
    return w;
}

static std::vector<std::string> split_tab(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    for (char c : line) {
        if (c == '\t') { out.push_back(cur); cur.clear(); }
        else if (c != '\r') cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

static void clean(std::string& s)   // single-line message for the protocol
{
    for (char& c : s) if (c == '\n' || c == '\r' || c == '\t') c = ' ';
}

int main(int argc, const char** argv)
{
    if (argc < 3) {
        std::cerr << "Usage: p1tool <sensor_profiles_dir> <jobs.tsv>\n";
        return 2;
    }
    P1::ImageSdk::SetSensorProfilesLocation(argv[1]);

    std::ifstream jobs(widen(argv[2]));
    if (!jobs) {
        std::cout << "@@ERR -1 cannot open job file" << std::endl;
        return 2;
    }
    std::string line;
    int idx = -1, nerr = 0;
    while (std::getline(jobs, line)) {
        if (line.empty() || line[0] == '#') continue;
        ++idx;
        std::vector<std::string> f = split_tab(line);
        if (f.size() < 12) {
            std::cout << "@@ERR " << idx << " bad job line (" << f.size() << " fields)" << std::endl;
            ++nerr;
            continue;
        }
        const std::string mode = f[0], iiqPath = f[1], outPath = f[2], cfgPath = f[3];
        try {
            const int x = std::stoi(f[4]), y = std::stoi(f[5]);
            const int w = std::stoi(f[6]), h = std::stoi(f[7]);
            const double scale = std::stod(f[8]);
            const int bits = std::stoi(f[9]), quality = std::stoi(f[10]);
            const bool lzw = std::stoi(f[11]) != 0;

            P1::ImageSdk::RawImage iiq(widen(iiqPath));          // read only
            P1::ImageSdk::ConvertConfig config = cfgPath.empty()
                ? P1::ImageSdk::ConvertConfig()
                : P1::ImageSdk::ConvertConfig::LoadFromFile(cfgPath.c_str());

            if (mode == "info") {
                P1::ImageSdk::Size cs = iiq.GetCanvasSize(config);
                std::cout << "@@INFO " << idx << " " << cs.width << " " << cs.height << " "
                          << iiq.Width() << " " << iiq.Height() << std::endl;
                continue;
            }
            const bool deep = bits == 16 && mode != "jpg";
            config.SetOutputFormat(deep ? P1::ImageSdk::BitmapFormat::Rgb48
                                        : P1::ImageSdk::BitmapFormat::Rgb24);
            // A settings file with "crop": true but no crop rectangle (e.g. APPN's
            // config_hr) develops an EMPTY bitmap for a whole frame, so a whole frame
            // is always requested as an explicit crop of the full canvas.
            if (w > 0 && h > 0) {
                config.SetCrop(x, y, w, h);
            } else {
                P1::ImageSdk::Size cs = iiq.GetCanvasSize(config);
                config.SetCrop(0, 0, cs.width, cs.height);
            }
            config.SetOutputScale(scale);

            P1::ImageSdk::BitmapImage bmp = iiq.Convert(config);
            const uint32_t bw = bmp.Width(), bh = bmp.Height();
            if (bw == 0 || bh == 0 || bmp.ByteSize() <= 0)
                throw std::runtime_error("the SDK returned an empty image");
            if (bmp.Format().GetPixelLayout() != P1::ImageSdk::PixelLayout::Interleaved)
                throw std::runtime_error("the SDK returned a non-interleaved image");

            if (mode == "tiff") {
                P1::ImageSdk::TiffConfig tc;
                if (lzw) tc.compression = P1::ImageSdk::TiffCompressionScheme_LZW;
                tc.commonConfig.includePreview = false;
                P1::ImageSdk::TiffWriter(widen(outPath), bmp, iiq, tc);
            } else if (mode == "jpg") {
                P1::ImageSdk::JpegConfig jc;
                jc.quality = quality;
                jc.commonConfig.includePreview = false;
                P1::ImageSdk::JpegWriter(outPath, bmp, iiq, jc);
            } else if (mode == "raw") {
                const size_t rowBytes = (size_t)bw * 3 * (deep ? 2 : 1);
                const size_t stride = bmp.Stride();
                std::shared_ptr<uint8_t> data = bmp.Data();
                FILE* fp = _wfopen(widen(outPath).c_str(), L"wb");
                if (!fp) throw std::runtime_error("cannot write " + outPath);
                for (uint32_t r = 0; r < bh; ++r)
                    fwrite(data.get() + (size_t)r * stride, 1, rowBytes, fp);
                fclose(fp);
            } else {
                throw std::runtime_error("unknown mode " + mode);
            }
            std::cout << "@@OK " << idx << " " << bw << " " << bh << " " << (deep ? 16 : 8) << std::endl;
        } catch (const std::exception& e) {
            std::string m = e.what();
            clean(m);
            std::cout << "@@ERR " << idx << " " << m << std::endl;
            ++nerr;
        }
    }
    return nerr ? 1 : 0;
}
