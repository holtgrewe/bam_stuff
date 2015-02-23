// ==========================================================================
//                            BAM Duplicate Masker
// ==========================================================================
// Copyright (c) 2015, Manuel Holtgrewe, Charite Universitaetsmedizin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@charite.de>
// ==========================================================================

#ifndef BAM_MARKDUP_SRC_BAM_MARKDUP_OPTIONS_H_
#define BAM_MARKDUP_SRC_BAM_MARKDUP_OPTIONS_H_

#include <iosfwd>
#include <string>
#include <stdexcept>
#include <vector>

// ----------------------------------------------------------------------------
// Class InvalidArgumentsException
// ----------------------------------------------------------------------------

class InvalidCommandLineArgumentsException : public std::runtime_error
{
public:
    InvalidCommandLineArgumentsException() : std::runtime_error("")
    {}
};

// ----------------------------------------------------------------------------
// Class BamMarkDupOptions
// ----------------------------------------------------------------------------

class BamMarkDupOptions
{
public:
    // Verbosity: 0 - quiet, 1 - normal, 2 - verbose, 3 - very verbose.
    int verbosity;

	// Command line as coming from argv.
	std::vector<std::string> commandLine;

    // Input alignment file (.bam).
    std::string inPath;
    // Output alignment file.
    std::string outPath;

    // Maximal expected template size (ignore entries with larger TLEN).
    int windowLength;
    // Force treating of paired reads as single-end reads.
    bool treatPairedAsSingle;

    BamMarkDupOptions() : verbosity(1), windowLength(4000), treatPairedAsSingle(false)
    {}

    void print(std::ostream & out) const;
};

// ----------------------------------------------------------------------------
// Function parseCommandLine()
// ----------------------------------------------------------------------------

BamMarkDupOptions parseCommandLine(int argc, char ** argv);

#endif  // #ifndef BAM_MARKDUP_SRC_BAM_MARKDUP_OPTIONS_H_
