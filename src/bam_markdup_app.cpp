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

#include "bam_markdup_app.h"

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include <seqan/bam_io.h>
#include <seqan/seq_io.h>
#include <seqan/simple_intervals_io.h>

#include "bam_markdup_options.h"
#include "progress_indicator.h"
#include "version.h"

namespace {  // anonymous namespace

// ---------------------------------------------------------------------------
// Class MyHash
// ---------------------------------------------------------------------------

struct MyHash
{
public:
    template <typename T, typename U>
    std::size_t operator()(std::pair<T, U> const & x) const
    {
        return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
    }

    std::size_t hashCode(seqan::CharString const & str) const
    {
        std::size_t hash = 0, strlen = length(str), i;
        if (strlen == 0)
            return hash;
        for (i = 0; i < strlen; i++)
            hash = (31 * hash) + str[i];
        return hash;
    }

    template <typename U, typename V>
    std::size_t operator()(std::tuple<seqan::CharString, U, V> const & x) const
    {
        return hashCode(std::get<0>(x)) ^ std::hash<U>()(std::get<1>(x)) ^ std::hash<V>()(std::get<2>(x));
    }

    template <typename T, typename U, typename V, typename W>
    std::size_t operator()(std::tuple<T, U, V, W> const & x) const
    {
        return std::hash<T>()(std::get<0>(x)) ^ std::hash<T>()(std::get<1>(x)) ^ std::hash<U>()(std::get<2>(x)) ^ std::hash<V>()(std::get<3>(x));
    }

};

// ---------------------------------------------------------------------------
// Function fileSize()
// ---------------------------------------------------------------------------

__uint64 fileSize(char const * path)
{
    std::ifstream in(path, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

// ---------------------------------------------------------------------------
// Function alignmentsStartInWindow()
// ---------------------------------------------------------------------------

// return true if the two alignments start within a window of the given size
bool alignmentsStartInWindow(seqan::BamAlignmentRecord const & lhs,
                             seqan::BamAlignmentRecord const & rhs,
                             int windowSize)
{
    if (lhs.rID != rhs.rID)
        return false;
    if (lhs.rID == seqan::BamAlignmentRecord::INVALID_REFID ||
        rhs.rID == seqan::BamAlignmentRecord::INVALID_REFID)
        return false;
    return (abs(lhs.beginPos - rhs.beginPos) < windowSize);
}

// ---------------------------------------------------------------------------
// Function isLeftOf()
// ---------------------------------------------------------------------------

// Original comparison function for strings by Heng Li from samtools.

inline int strnum_cmp(const char *a, const char *b)
{
    char *pa, *pb;
    pa = (char*)a; pb = (char*)b;
    while (*pa && *pb) {
        if (isdigit(*pa) && isdigit(*pb)) {
            long ai, bi;
            ai = strtol(pa, &pa, 10);
            bi = strtol(pb, &pb, 10);
            if (ai != bi) return ai<bi? -1 : ai>bi? 1 : 0;
        } else {
            if (*pa != *pb) break;
            ++pa; ++pb;
        }
    }
    if (*pa == *pb)
        return (pa-a) < (pb-b)? -1 : (pa-a) > (pb-b)? 1 : 0;
    return *pa<*pb? -1 : *pa>*pb? 1 : 0;
}

// with tie-breaking by name
bool isLeftOf(seqan::BamAlignmentRecord const & lhs, seqan::BamAlignmentRecord const & rhs)
{
    return ((unsigned)lhs.rID < (unsigned)rhs.rID ||
            (lhs.rID == rhs.rID && lhs.beginPos < rhs.beginPos) ||
            (lhs.rID == rhs.rID && lhs.beginPos == rhs.beginPos && strnum_cmp(toCString(lhs.qName), toCString(rhs.qName)) < 0));
}

// without tie-breaking by name
bool isLeftOfPosOnly(seqan::BamAlignmentRecord const & lhs, seqan::BamAlignmentRecord const & rhs)
{
    return ((unsigned)lhs.rID < (unsigned)rhs.rID || (lhs.rID == rhs.rID && lhs.beginPos < rhs.beginPos));
}

// ---------------------------------------------------------------------------
// Class DuplicateBamProcessor
// ---------------------------------------------------------------------------

/*!
 * Stream through BamAlignmentRecord objects and execute callback on batches of those aligning on duplicate positions.
 */

class DuplicateBamProcessor
{
public:
    DuplicateBamProcessor(seqan::BamFileOut & bamFileOut, BamMarkDupOptions const & options) :
        bamFileOut(bamFileOut), options(options)
    {}

    ~DuplicateBamProcessor();

    void push(std::vector<seqan::BamAlignmentRecord *> & records);

private:

    // process buffer; if isFinal then also write out the last segment
    void process(bool isFinal = false);

    typedef std::vector<seqan::BamAlignmentRecord *>::iterator TIter;
    typedef std::pair<TIter, TIter> TRange;

    // processing of single-ended reads
    void processSE(TRange window);
    // processing of paired-ended reads
    void processPE(TRange window);

    // partition the given range of records at the same start position and
    // call markSingleEndDupes on each batch
    void partitionAndMarkDupesSE(std::vector<seqan::BamAlignmentRecord *>::iterator itBegin,
                                 std::vector<seqan::BamAlignmentRecord *>::iterator itEnd);
    // mark single-ended duplicates in already selected batch of duplicates
    void markSingleEndDupes(std::vector<seqan::BamAlignmentRecord *> & duplicates);
    // mark paired-end duplicates in already selected batch of duplicates
    void markPairedEndDupes(
            std::vector<std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *>> & duplicates);

    // buffer with yet-to-be-processed records
    std::vector<seqan::BamAlignmentRecord *> buffer;

    // sink for processed BAM files
    seqan::BamFileOut & bamFileOut;
    // configuration
    BamMarkDupOptions options;
};

void DuplicateBamProcessor::partitionAndMarkDupesSE(std::vector<seqan::BamAlignmentRecord *>::iterator itBegin,
                                                    std::vector<seqan::BamAlignmentRecord *>::iterator itEnd)
{
    // partition the given range by (begin pos, end pos)
    std::unordered_map<std::pair<int, int>, std::vector<seqan::BamAlignmentRecord *>, MyHash> partitions;
    for (auto it = itBegin; it != itEnd; ++it)
    {
        if (hasFlagMultiple(**it) && !options.treatPairedAsSingle)
            continue;  // ignore single-end reads
        if (hasFlagSupplementary(**it))
            continue;  // ignore supplementary alignments
        int const beginPos = (*it)->beginPos;
        int const endPos = beginPos + getAlignmentLengthInRef(**it);
        partitions[std::make_pair(beginPos, endPos)].push_back(*it);
    }

    // mark duplicates in each partition
    for (auto & entry : partitions)
        markSingleEndDupes(entry.second);
}

void DuplicateBamProcessor::markSingleEndDupes(std::vector<seqan::BamAlignmentRecord *> & duplicates)
{
    if (duplicates.size() <= 1u)
        return;

    // obtain best by quality
    auto cmpMapQ = [](seqan::BamAlignmentRecord const * lhs, seqan::BamAlignmentRecord const * rhs) {
        return lhs->mapQ < rhs->mapQ;
    };
    auto best = std::max_element(duplicates.begin(), duplicates.end(), cmpMapQ);
    if (options.verbosity >= 3)
        std::cerr << "best read: " << (*best)->qName << "/" << hasFlagLast(**best) + 1 << "\n";

    // mask all but best with appropriate BAM flag
    for (auto it = duplicates.begin(); it != duplicates.end(); ++it)
    {
        if (options.verbosity >= 3)
            std::cerr << "marking as SE duplicate: " << (*it)->qName << "/" << hasFlagLast(**it) + 1 << "\n";
        if (it != best)
            (*it)->flag |= seqan::BAM_FLAG_DUPLICATE;
    }
}

void DuplicateBamProcessor::markPairedEndDupes(
        std::vector<std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *>> & duplicates)
{
    if (duplicates.size() <= 1u)
        return;

    // obtain best by quality
    auto cmpMapQ = [](std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *> const & lhs,
                      std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *> const & rhs) {
        return (lhs.first->mapQ + lhs.second->mapQ) < (rhs.first->mapQ + rhs.second->mapQ);
    };
    auto best = std::max_element(duplicates.begin(), duplicates.end(), cmpMapQ);
    if (options.verbosity >= 3)
        std::cerr << "best pair " << best->first->qName << "\n";

    // mask all but best with appropriate BAM flag
    for (auto it = duplicates.begin(); it != duplicates.end(); ++it)
        if (it != best) {
            if (options.verbosity >= 3)
                std::cerr << "marking as PE duplicate: " << it->first->qName << "/" << hasFlagLast(*it->second) + 1 << "\n";
            it->first->flag |= seqan::BAM_FLAG_DUPLICATE;
            it->second->flag |= seqan::BAM_FLAG_DUPLICATE;
        }
}

DuplicateBamProcessor::~DuplicateBamProcessor()
{
    // process buffer, setting isFinal flag to true
    process(true);
}

void DuplicateBamProcessor::push(std::vector<seqan::BamAlignmentRecord *> & records)
{
    // move records to buffer
    buffer.insert(buffer.end(), records.begin(), records.end());
    records.clear();
    // process buffer
    process();
}

void DuplicateBamProcessor::process(bool isFinal)
{
    if (buffer.empty())
        return;  // guard against empty buffer
    if (options.verbosity >= 2)
        std::cerr << "buffer size == " << buffer.size() << "\n";

    // Get longest prefix of buffer, such that the rightmost read is more than options.windowLength from the last
    // read in buffer. If this is the final alignment round (and no more alignments are expected) then consider
    // all of buffer.
    auto range = std::make_pair(buffer.begin(), buffer.begin());
    if (isFinal)
        range.second = buffer.end();
    else
        while (range.second != buffer.end() &&
                !alignmentsStartInWindow(**range.second, **(buffer.end() - 1), options.windowLength))
            ++range.second;

    // Process single-end and paired-end records in range (options.windowLength > 1 ensures that SE processing
    // works correctly).
    processSE(range);
    if (!options.treatPairedAsSingle)
        processPE(range);

    // write out all processed, free them, and remove them from buffer
    for (auto it = range.first; it != range.second; ++it)
    {
        writeRecord(bamFileOut, **it);
        delete *it;
    }
    buffer.erase(range.first, range.second);
}

void DuplicateBamProcessor::processSE(TRange window)
{
    auto itBegin = window.first;

    while (itBegin != window.second)
    {
        if (hasFlagUnmapped(**itBegin))
            break;  // reached first unaligned record

        // get range of records at equal start positions
        auto range = std::make_pair(itBegin, itBegin);
        while (range.second != window.second && !isLeftOfPosOnly(**range.first, **range.second))
            ++range.second;
        if (range.first == range.second)
            break;  // range was empty
        if (!isLeftOfPosOnly(**range.first, **(range.second - 1)))
            break;  // last not greater than first

        partitionAndMarkDupesSE(range.first, range.second);
        itBegin = range.second;
    }
}

void DuplicateBamProcessor::processPE(TRange window)
{
    if (options.verbosity >= 3)
        std::cerr << "PE WINDOW\n\n";

    // obtain records, sorted by (qName, beginPos, firstFlag)
    std::unordered_map<std::tuple<seqan::CharString, bool, int>, seqan::BamAlignmentRecord *, MyHash> records;
    for (auto elem : buffer)
        if (hasFlagMultiple(*elem)
                && !hasFlagSupplementary(*elem)
                && (!hasFlagUnmapped(*elem) || !hasFlagNextUnmapped(*elem))
                && (elem->rID == elem->rNextId))
        {
            auto key = std::make_tuple(elem->qName, hasFlagLast(*elem), elem->beginPos);
            if (records.find(key) != records.end())
                std::cerr << "WARNING: found duplicate alignment for the same read " << elem->qName
                          << "/" << (hasFlagLast(*elem)) << "\n";
            records[key] = elem;
            if (options.verbosity >= 3)
                std::cerr << "(" << std::get<0>(key) << ", " << std::get<1>(key) << ", " << std::get<2>(key)
                          << ") -> " << elem->qName << "/" << (hasFlagLast(*elem) + 1) << "\n";
        }

    // collect records, obtaining pairing
    std::vector<std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *>> pairs;
    for (auto it = window.first; it != window.second; ++it)
        if (hasFlagMultiple(**it)  // paired
                && !hasFlagSupplementary(**it)
                && (!hasFlagUnmapped(**it) || !hasFlagNextUnmapped(**it))  // one mate aligned
                && ((*it)->rID == (*it)->rNextId)  // aligned on same ref
                && (std::make_pair((*it)->beginPos, hasFlagLast(**it))  // left, ties broken by first flag
                    < std::make_pair((*it)->pNext, hasFlagLast(**it))))
        {
            auto key = std::make_tuple((*it)->qName, !hasFlagLast(**it), (*it)->pNext);
            if (options.verbosity >= 3)
                std::cerr << "Looking for (" << std::get<0>(key) << ", " << std::get<1>(key) << ", "
                          << std::get<2>(key) << ")\n";
            auto other = records.find(key);
            if (other == records.end())
            {
                if (abs((*it)->pNext - (*it)->beginPos) < options.windowLength)
                    std::cerr << "WARNING: could not find mate for " << (*it)->qName << "/" << (hasFlagLast(**it) + 1)
                              << " (broken BAM file?)\n";
                else if (options.verbosity >= 3)
                    std::cerr << "INFO: ignoring "  << (*it)->qName << "/" << (hasFlagLast(**it) + 1)
                              << " since mate is too far away.\n";
            }
            else if ((*it)->beginPos < other->second->beginPos)
            {
                pairs.push_back(std::make_pair(*it, other->second));
            }
            else
            {
                pairs.push_back(std::make_pair(other->second, *it));
            }
        }

    // partition buffer, adding the necessary flags
    typedef std::tuple<int, int, int, int> TKey;
    typedef std::pair<seqan::BamAlignmentRecord *, seqan::BamAlignmentRecord *> TValue;
    std::unordered_map<TKey, std::vector<TValue>, MyHash> partitions;
    for (auto pair : pairs)
    {
        TKey key(pair.first->beginPos, pair.first->beginPos + getAlignmentLengthInRef(*pair.first),
                 pair.second->beginPos, pair.second->beginPos + getAlignmentLengthInRef(*pair.second));
        partitions[key].push_back(pair);
    }
    // apply marks
    for (auto & elem : partitions)
        markPairedEndDupes(elem.second);
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// Class BamMarkDupAppImpl
// ---------------------------------------------------------------------------

class BamMarkDupAppImpl
{
public:
    BamMarkDupAppImpl(BamMarkDupOptions const & options) :
        options(options), bamFileIn(bamFileOut)
    {}

    void run();

private:

    // Process regions one-by-one.
    void processRecords();

    // Open files.
    void openFiles();

    // Program configuration.
    BamMarkDupOptions options;

    // Objects used for I/O.
    seqan::BamFileOut bamFileOut;
    seqan::BamFileIn bamFileIn;

    // BAM header is read into this variable.
    seqan::BamHeader bamHeader;
};

void BamMarkDupAppImpl::run()
{
    double startTime = seqan::sysTime();

    if (options.verbosity >= 1)
    {
        std::cerr << "BAM Duplicate Masker\n"
                  << "====================\n\n";
        options.print(std::cerr);
    }

    // Open Files
    openFiles();

    // Process BAM Records

    if (options.verbosity >= 1)
        std::cerr << "\n"
                  << "__PROCESSING RECORDS_____________________________________________\n"
                  << "\n";

    processRecords();

    if (options.verbosity >= 1)
        std::cerr << "\nTotal time: " << (seqan::sysTime() - startTime) << " s\n"
                  << "All done. Have a nice day!\n";
}

void BamMarkDupAppImpl::openFiles()
{
    if (options.verbosity >= 1)
    {
        std::cerr << "\n"
                  << "__OPENING FILES__________________________________________________\n"
                  << "\n";
    }

    std::cerr << "Opening " << options.inPath << " ...";
    if (options.inPath == "-" && !open(bamFileIn, std::cin))
        throw new std::runtime_error("Could not stdin for input.");
    else if (options.inPath != "-" && !open(bamFileIn, options.inPath.c_str()))
        throw new std::runtime_error("Could not open input file.");
    std::cerr << " OK\n";

    readHeader(bamHeader, bamFileIn);
    seqan::BamHeaderRecord pgRecord;
    pgRecord.type = seqan::BAM_HEADER_PROGRAM;
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("ID", "bam_markdup"));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("PN", "bam_markdup"));
    std::string commandLine;
    for (auto token : options.commandLine)
    {
        if (!commandLine.empty())
            commandLine.append(" ");
        commandLine.append(token);
    }
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("VN", PROGRAM_VERSION));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("CL", commandLine.c_str()));
    appendValue(pgRecord.tags, seqan::BamHeaderRecord::TTag("DS", "masking of duplicate alignments"));
    appendValue(bamHeader, pgRecord);

    std::cerr << "Opening " << options.outPath << " ...";
    if (!open(bamFileOut, options.outPath.c_str()))
        throw new std::runtime_error("Could not open output file.");
    std::cerr << " OK\n";
    writeHeader(bamFileOut, bamHeader);
}

void BamMarkDupAppImpl::processRecords()
{
    unsigned const MIB = 1; //1024 * 1024;
    double startTime = seqan::sysTime();

    std::unique_ptr<ProgressBar> progress;
    progress.reset(new ProgressBar(std::cerr, 0, fileSize(options.inPath.c_str()) / MIB,
                                   (options.verbosity == 1)));
    progress->setLabel(options.inPath.c_str());
    progress->updateDisplay();

    DuplicateBamProcessor processor(bamFileOut, options);
    std::vector<seqan::BamAlignmentRecord *> chunk;
    size_t const CHUNK_SIZE = 10000;
    while (!atEnd(bamFileIn))
    {
        // read in chunk and put into processor
        while (chunk.size() < CHUNK_SIZE && !atEnd(bamFileIn))
        {
            std::unique_ptr<seqan::BamAlignmentRecord> ptr(new seqan::BamAlignmentRecord);
            readRecord(*ptr, bamFileIn);
            chunk.push_back(ptr.release());
        }
        if (options.verbosity >= 2)
            std::cerr << "reader at " << chunk.front()->rID << ":" << chunk.front()->beginPos << " (chunk size == " << chunk.size() << ")\n";
        processor.push(chunk);

        // update progress bar
        progress->advanceTo((position(bamFileIn) >> 16) / MIB);
    }

    progress->finish();

    if (options.verbosity >= 1)
        std::cerr << "Masking time: " << (seqan::sysTime() - startTime) << " s\n";
}

// ---------------------------------------------------------------------------
// Class BamMarkDupApp
// ---------------------------------------------------------------------------

BamMarkDupApp::BamMarkDupApp(BamMarkDupOptions const & options) :
        impl(new BamMarkDupAppImpl(options))
{}

BamMarkDupApp::~BamMarkDupApp()  // for pimpl
{}

void BamMarkDupApp::run()
{
    impl->run();
}
