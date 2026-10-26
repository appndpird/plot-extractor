# Third-party notices

## APPN-UQ PhaseOne_Pipeline

PlotExtractor's IIQ plot extraction follows the approach of the APPN-UQ PhaseOne_Pipeline
(Australian Plant Phenomics Network, The University of Queensland node; contact
appn@uq.edu.au). `iiq_tools/src/p1tool.cpp` is adapted from its `IIQ_Proc` tools
(`IIQ2TIFF.cpp`, `IIQHRclip.cpp`), and `iiq_tools/presets/appn_hr.json` and
`iiq_tools/presets/appn_template.json` are its `config_hr` settings.

```
MIT License

Copyright (c) 2025 APPN-UQ

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

## Phase One ImageSDK

IIQ files are developed by the Phase One ImageSDK, Copyright Phase One A/S, used under
the Phase One SDK Software License Agreement. The SDK is not part of the PlotExtractor
source code. Where an application build includes the SDK runtime (`iiq_tools/sdk/`), it
is there under Phase One's licence (its `PhaseOne_SDK_LICENSE.md` and
`PhaseOne_SDK_CREDITS.md`, the third-party libraries used by the SDK, are included
alongside) and may not be extracted or used on its own.
