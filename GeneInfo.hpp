/*
 * =====================================================================================
 *
 *       Filename:  GeneInfo.hpp
 *
 *    Description:  A simple class storing and manipulating all the needed data from a 
 *                  gene annotation table.
 *
 *        Version:  1.0
 *        Created:  05/01/11 10:20:57
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jianxing Feng (), jianxing.tongji@gmail.com
 *        Company:  Tongji Univ.
 *
 * =====================================================================================
 */

#ifndef GeneInfo_H
#define GeneInfo_H

#include <vector>
#include <string>
#include <map>
#include <list>
#include <set>
#include <assert.h>
#include "Utility.hpp"

using namespace std;

class Interval
{
public:
    int mStart;
    int mEnd;
    void* mpData;     // useful to sort

    Interval(){}
    Interval(int start, int end, void* p_data = NULL){
        mStart = start; mEnd = end; mpData = p_data;
    }

    int Length(){return mEnd - mStart;}

    bool operator<(const Interval &other) const {
        return (mStart < other.mStart || mStart == other.mStart && mEnd < other.mEnd);
    }

    bool operator!=(const Interval &other) const {
        return (mStart != other.mStart || mEnd != other.mEnd || mpData != other.mpData);
    }
};

class IntervalTree
{
public:
    vector<Interval*> mIntervals;
    IntervalTree* mpLeft;
    IntervalTree* mpRight;
    int mCenter;
    
    IntervalTree(){mpLeft = mpRight = NULL;}

    void
    Print()
    {
        if (mpLeft) mpLeft->Print();
        for (size_t i = 0; i < mIntervals.size(); ++i)
            cout << mIntervals[i]->mStart << "\t" << mIntervals[i]->mStart << endl;
        if (mpRight) mpRight->Print();
    }

    void
    IntervalCount(int& count)
    {
        if (mpLeft) mpLeft->IntervalCount(count);
        count += mIntervals.size();
        if (mpRight) mpRight->IntervalCount(count);
    }

    void
    Init(const vector<Interval*>& intervals, size_t mincnt = 4)
    {
        mpLeft = mpRight = NULL;

        if (intervals.size() < mincnt)
        {
            mIntervals = intervals;
            return;
        }

        mCenter = (intervals[intervals.size()/2]->mStart + intervals[intervals.size()/2]->mEnd)/2;

        vector<Interval*> lefts, rights;
        for (size_t i = 0; i < intervals.size(); ++i)
        {
            if (intervals[i]->mEnd <= mCenter)
                lefts.push_back(intervals[i]);
            else if (intervals[i]->mStart > mCenter)
                rights.push_back(intervals[i]);
            else
                mIntervals.push_back(intervals[i]);
        }
        if (lefts.size() == intervals.size())
        {
            mIntervals = lefts;
            return;
        }
        if (rights.size() == intervals.size())
        {
            mIntervals = rights;
            return;
        }
        if (lefts.size() > 0)
        {
            mpLeft = new IntervalTree();
            mpLeft->Init(lefts, mincnt);
        }
        if (rights.size() > 0)
        {
            mpRight = new IntervalTree();
            mpRight->Init(rights, mincnt);
        }
    }

    ~IntervalTree()
    {
        if (mpLeft) delete mpLeft;
        if (mpRight) delete mpRight;
    }

    void Find(int pos, vector<Interval*>& overlapping)
    {
        for (size_t i = 0; i < mIntervals.size(); ++i)
            if (mIntervals[i]->mEnd > pos && mIntervals[i]->mStart <= pos)
                overlapping.push_back(mIntervals[i]);

        if (mpLeft && pos < mCenter)
            mpLeft->Find(pos, overlapping);
        if (mpRight && pos > mCenter)
            mpRight->Find(pos, overlapping);
    }

    void Find(Interval& inter, vector<Interval*>& overlapping)
    {
        for (size_t i = 0; i < mIntervals.size(); ++i)
            if (mIntervals[i]->mEnd > inter.mStart && mIntervals[i]->mStart < inter.mEnd) 
                overlapping.push_back(mIntervals[i]);

        if (mpLeft && inter.mStart < mCenter)
            mpLeft->Find(inter, overlapping);
        if (mpRight && inter.mEnd > mCenter)
            mpRight->Find(inter, overlapping);
    }
};

class Gene: public Interval
{
public:
    string mChr;
    char mStrand;
    string mGeneSymbol;
    string mTranscriptID;
    vector<Interval> mExons;

    /* All the exons overlappaed are merged */
    void UnionExons()
    {
        sort(mExons.begin(), mExons.end());
     
        int idx = 0;
        for (size_t i = 1; i < mExons.size(); ++i)
        {
            if (mExons[i].mStart < mExons[idx].mEnd)
            {
                if (mExons[i].mEnd > mExons[idx].mEnd)     // In case an interval is completely included in another one
                    mExons[idx].mEnd = mExons[i].mEnd;
            }
            else
                mExons[++idx] = mExons[i];
        }
        mExons.resize(idx+1);
     
        mStart = mExons[0].mStart;
        mEnd = mExons[mExons.size()-1].mEnd;
    }
};

// For each segment/junction, this structure stores all the information including mapped reads and refsequence.
class ExtraInfo
{
public:
    ExtraInfo(){mCount = 0;}

    list<Gene*> mGenes;
    list<int> mMapPos;
    list<int> mMapFlag; 
    int mSecondPart;              // Used only for junction
    string mFirstPartSeq;         // Used only for junction
    int mCount;
};

class GeneInfo
{
typedef map<string, Gene> map_str2gene_t;
typedef map<string, list<Interval> > map_str2list_inter_t;

public:
    map_str2gene_t mAllGenes;
    map_str2gene_t mAllTranscripts;
    map_str2list_inter_t mAllJunctions;
    map_str2list_inter_t mAllSegments;
    map<string, IntervalTree> mSegmentTree;
    map<string, IntervalTree> mJunctionTree;
    int mFlankLen;       // The length of flank sequences for each junction. The maximum length should be the read length
    bool mbStrandSpecific;
    bool mbCountJunc;    // In case only gene read count is needed, junction tree is no need to be scanned when loading short reads.

    // For loading refseq
    string mLastChr;
    int mLastPos;
    vector<Interval> mAllBedList;
    int mLastBedIdx;
    list<Interval> mTailList;
    int mSegCnt;
    int mJuncCnt;
    int mRefSeqBufLen;
    string mRefSeqBuf;
    int mLongestLine;

    // For output
    int mVerbosLevel;
    string mOutputFile;
    ofstream* mpOutput;
    
    int mSenseReadCount;
    int mAntiSenseReadCount;

    GeneInfo(string output_rawinfo, int flank_len, bool b_strand_specific, bool b_count_junc = true, int verbos_level = 2)
    {
        mOutputFile = output_rawinfo;
        mFlankLen = flank_len;
        mbStrandSpecific = b_strand_specific;
        mbCountJunc = b_count_junc;
        mVerbosLevel = verbos_level;
    }

    ~GeneInfo()
    {
        // After FinishLoadingExons is called. mAllSegments and mAllJunctions 
        // contains created objects
         for_each_ele_in_group(iter, map_str2list_inter_t, mAllSegments)
        {
            const string& chr = iter->first;
            list<Interval>& curr_segs = iter->second;
            for_each_ele_in_group(sub_iter, list<Interval>, curr_segs)
            {
                ExtraInfo* p_extra = (ExtraInfo*)(sub_iter->mpData);
                delete p_extra;
            }
            list<Interval>& curr_juncs = mAllJunctions[chr];
            for_each_ele_in_group(sub_iter, list<Interval>, curr_juncs)
            {
                ExtraInfo* p_extra = (ExtraInfo*)(sub_iter->mpData);
                delete p_extra;
            }
        }
    }

    void OnAnExon(const string& chr, char strand, const string& gene_id, const string& tran_id, int start, int end)
    {
        if (mAllTranscripts.find(tran_id) == mAllTranscripts.end())
        {
            Gene a_gene;
            a_gene.mChr = chr; 
            if (mbStrandSpecific)
                a_gene.mChr += strand;
            a_gene.mStrand = strand;
            a_gene.mGeneSymbol = gene_id;
            a_gene.mTranscriptID = tran_id;
            mAllTranscripts[tran_id] = a_gene;
        }

        Gene& a_gene = mAllTranscripts[tran_id];
        a_gene.mExons.push_back(Interval(start, end));
    } // OnAnExon

    void FinishLoadingExons()
    {
        typedef map<string, vector<int> > map_str2vec_int_t;
        typedef map<string, set<int> > map_str2set_int_t;
        map_str2set_int_t chr2breaks_set;
        // It is possible that a gene (identified by gene symbol) has several transcripts.
        // It is possible that a gene (identified by gene symbol) has different transcripts on different strand / chromosome.
        // For a junction, it should only be transcript specific.
        // If the data is not strand specific, a junction or segment with the same positions should be treated as duplicates.

        // Generate mAllGenes based on mAllTranscripts
        for_each_ele_in_group(iter, map_str2gene_t, mAllTranscripts)
        {
            Gene& a_gene = iter->second;
            if (mAllGenes.find(a_gene.mGeneSymbol) == mAllGenes.end())
            {
                Gene new_gene;
                new_gene.mChr = a_gene.mChr;
                new_gene.mStrand = a_gene.mStrand;
                new_gene.mGeneSymbol = a_gene.mGeneSymbol;
                new_gene.mTranscriptID = a_gene.mTranscriptID;
                mAllGenes[a_gene.mGeneSymbol] = new_gene;
            }

            Gene& new_gene = mAllGenes[a_gene.mGeneSymbol];
            // If two transcripts of the same gene appear at different chromosome or different strand,
            // keep only the first transcript
            if (new_gene.mChr == a_gene.mChr && new_gene.mStrand == a_gene.mStrand)
                new_gene.mExons.insert(new_gene.mExons.end(), a_gene.mExons.begin(), a_gene.mExons.end());
        }
        
        // Get all the junctions based on mAllTranscripts and being aware of strand specificity
        for_each_ele_in_group(iter, map_str2gene_t, mAllTranscripts)
        {
            string chr = iter->second.mChr;
            if (mAllJunctions.find(chr) == mAllJunctions.end())
            {
                mAllJunctions[chr] = list<Interval>();
                chr2breaks_set[chr] = set<int>();
            }
            list<Interval>& curr_juncs = mAllJunctions[chr];
            set<int>& curr_breaks = chr2breaks_set[chr];

            vector<Interval>& exons = iter->second.mExons;
            sort(exons.begin(), exons.end());
            for (size_t i = 1; i < exons.size(); ++i)
                // Note that a pointer to gene instead of transcript is used here, because we care about
                // expression levels of genes instead of transcripts.
                curr_juncs.push_back(Interval(exons[i-1].mEnd, exons[i].mStart, &mAllGenes[iter->second.mGeneSymbol]));

            for (size_t i = 0; i < exons.size(); ++i)
            {
                curr_breaks.insert(exons[i].mStart); 
                curr_breaks.insert(exons[i].mEnd); 
            }
        }
        
        // Remove duplicated junctions with the same start, end positions and mpData
        // Put all elements with the same start and end positions but different mpData
        // into a vector
        for_each_ele_in_group(iter, map_str2list_inter_t, mAllJunctions)
        {
            list<Interval>& curr_juncs = iter->second;
            curr_juncs.sort();

            list<Interval>::iterator new_iter = curr_juncs.begin();
            for_each_ele_in_group(junc_iter, list<Interval>, curr_juncs)
                if (*new_iter != *junc_iter)
                    *(++new_iter) = *junc_iter;
            curr_juncs.erase(++new_iter, curr_juncs.end());

            new_iter = curr_juncs.begin();
            for_each_ele_in_group(junc_iter, list<Interval>, curr_juncs)
            {
                Gene* p_curr_gene = (Gene*)(junc_iter->mpData);

                bool b_new_data = false;
                if (junc_iter == curr_juncs.begin())
                {
                    new_iter->mpData = new ExtraInfo;
                    b_new_data = true;
                }
                if (new_iter->mStart != junc_iter->mStart || new_iter->mEnd != junc_iter->mEnd)
                {
                    (++new_iter)->mpData = new ExtraInfo;
                    b_new_data = true;
                }

                ExtraInfo& extra = *((ExtraInfo*)(new_iter->mpData));
                extra.mGenes.push_back(p_curr_gene);
                if (b_new_data)
                {
                    new_iter->mStart = junc_iter->mStart; 
                    new_iter->mEnd = junc_iter->mEnd; 
                }
            }
            curr_juncs.erase(++new_iter, curr_juncs.end());

            // Adjust positions to be ready for scanning short reads
            for_each_ele_in_group(junc_iter, list<Interval>, curr_juncs)
            {
                ExtraInfo& extra = *((ExtraInfo*)(junc_iter->mpData));
                extra.mSecondPart = junc_iter->mEnd;                // TRICKY
                junc_iter->mEnd = junc_iter->mStart;                // TRICKY
                junc_iter->mStart = junc_iter->mEnd - mFlankLen;    // TRICKY
                if (junc_iter->mStart < 0)  
                    junc_iter->mStart = 0;  
            }
        }

        // Build interval trees for junctions 
        for_each_ele_in_group(iter, map_str2list_inter_t, mAllJunctions)
        {
            const string& chr = iter->first;
            list<Interval>& curr_juncs_list = mAllJunctions[chr];
            vector<Interval*> curr_juncs_vec;
            curr_juncs_vec.reserve(curr_juncs_list.size());
            for_each_ele_in_group(junc_iter, list<Interval>, curr_juncs_list)
                curr_juncs_vec.push_back(&(*junc_iter));
            mJunctionTree[chr] = IntervalTree();
            mJunctionTree[chr].Init(curr_juncs_vec);
        }

        // Processing segments...
        // Get all the collapsed exons and group them by chromosome
        for_each_ele_in_group(iter, map_str2gene_t, mAllGenes)
            iter->second.UnionExons();

        // Assuming that all the transcripts of a gene appear at the same strand
        // of some chromosome.
        for_each_ele_in_group(iter, map_str2gene_t, mAllGenes)
        {
            string chr = iter->second.mChr;
            if (mAllSegments.find(chr) == mAllSegments.end())
                mAllSegments[chr] = list<Interval>();
            list<Interval>& curr_segs = mAllSegments[chr];
            for_each_ele_in_group(seg_iter, vector<Interval>, (iter->second).mExons)
                curr_segs.push_back(Interval(seg_iter->mStart, seg_iter->mEnd, &(iter->second)));
        }

        int all_seg_cnt = 0;
        // Break all the exons into segments
        for_each_ele_in_group(iter, map_str2set_int_t, chr2breaks_set)
        {
            const string& chr = iter->first;
            list<Interval> curr_segs = mAllSegments[chr];     // Make a copy
            vector<Interval> curr_breaks;
            vector<Interval*> curr_breaks_pointer;
            curr_breaks.resize(iter->second.size());
            curr_breaks_pointer.resize(iter->second.size());
            int idx = 0;
            for_each_ele_in_group(break_iter, set<int>, iter->second)
            {
                curr_breaks[idx].mStart = *break_iter;
                curr_breaks[idx].mEnd = curr_breaks[idx].mStart + 1;
                curr_breaks_pointer[idx] = &curr_breaks[idx];
                ++idx;
            }
            IntervalTree break_tree;
            break_tree.Init(curr_breaks_pointer);
            
            list<Interval>& new_segs = mAllSegments[chr];
            new_segs.resize(0);
            for_each_ele_in_group(seg_iter, list<Interval>, curr_segs)
            {
                vector<Interval*> overlapping;
                break_tree.Find(*seg_iter, overlapping);
                sort(overlapping.begin(), overlapping.end());

                int start = seg_iter->mStart;
                int end = seg_iter->mEnd;

                assert(overlapping.size() > 0);
                assert(start == overlapping[0]->mStart);
                
                for (size_t i = 1; i < overlapping.size(); ++i)
                {
                    new_segs.push_back(Interval(start, overlapping[i]->mStart, seg_iter->mpData));
                    start = overlapping[i]->mStart;
                }
                new_segs.push_back(Interval(start, end, seg_iter->mpData));
            }
            //cout << chr << "\t" << curr_segs.size() << " -> " << new_segs.size();

            // Put gene information of segments with the same start and end positions but different mpData
            // into a vector
            new_segs.sort();
            list<Interval>::iterator new_iter = new_segs.begin();
            for_each_ele_in_group(seg_iter, list<Interval>, new_segs)
            {
                Gene* p_curr_gene = (Gene*)(seg_iter->mpData);
                bool b_new_data = false;
                if (seg_iter == new_segs.begin())
                {
                    new_iter->mpData = new ExtraInfo;
                    b_new_data = true;
                }
                if (new_iter->mStart != seg_iter->mStart || new_iter->mEnd != seg_iter->mEnd)
                {
                    (++new_iter)->mpData = new ExtraInfo;
                    b_new_data = true;
                }
                if (b_new_data)
                {
                    new_iter->mStart = seg_iter->mStart; 
                    new_iter->mEnd = seg_iter->mEnd; 
                }
                //cout << &(*new_iter) << "\t" << &(*seg_iter) << "\t" 
                //     << new_iter->mStart << "\t" << new_iter->mEnd << "\t" << seg_iter->mStart << "\t" << seg_iter->mEnd << endl;
                ExtraInfo& extra = *((ExtraInfo*)(new_iter->mpData));
                extra.mGenes.push_back(p_curr_gene);
            }
            new_segs.erase(++new_iter, new_segs.end());

            //cout << " -> " << new_segs.size() << endl;
            all_seg_cnt += new_segs.size();
        }
        //cout << "all_seg_cnt = " << all_seg_cnt << endl;

        // Build interval trees for segments
        for_each_ele_in_group(iter, map_str2list_inter_t, mAllSegments)
        {
            const string& chr = iter->first;
            list<Interval>& curr_segs = iter->second;
            vector<Interval*> curr_segs_pointer;
            curr_segs_pointer.reserve(curr_segs.size());
            for_each_ele_in_group(seg_iter, list<Interval>, curr_segs)
                curr_segs_pointer.push_back(&(*seg_iter));
            mSegmentTree[chr] = IntervalTree();
            mSegmentTree[chr].Init(curr_segs_pointer);
        }
    } // FinishLoadingExons

    void BeginLoadingShortReads()
    {
        mSenseReadCount = 0;
        mAntiSenseReadCount = 0;
    }

    /* Call this function after gene annotation has been loaded */
    // Assuming that start position is 0-based, which is different from SAM format
    void OnAShortRead(string chr, int start, int flag, int read_len, const string& match_string)
    {
        bool b_positive = !(flag & 16);
        bool b_junc = (match_string.find('N') != string::npos);

        char strand = '-';
        if (b_positive)
            strand = '+';
        if (mbStrandSpecific)
            chr += strand;

        if (b_junc && mbCountJunc)
        {
            IntervalTree& tree = mJunctionTree[chr];
            vector<Interval*> overlapping;
            tree.Find(start, overlapping);

            size_t M_pos = match_string.find_first_of('M');
            size_t N_pos = match_string.find_first_of('N');
            string ts = match_string.substr(M_pos + 1, N_pos - M_pos - 1);
            int gap_len = atoi(ts.data());
            for (size_t i = 0; i < overlapping.size(); ++i)
            {
                ExtraInfo& extra = *((ExtraInfo*)(overlapping[i]->mpData));
                if (overlapping[i]->mEnd + gap_len == extra.mSecondPart)
                {
                    if (start - overlapping[i]->mStart >= 2 * mFlankLen)
                        cerr << "WARNING: Flank length for a junction is too short, please set a value at least " << (start - overlapping[i]->mStart + 1) / 2 << endl;
                    extra.mMapPos.push_back(start - overlapping[i]->mStart);   
                    extra.mMapFlag.push_back(flag);
                }
            } 
        }
        else
        {
            if (!b_positive) 
                start += read_len - 1;

            IntervalTree& tree = mSegmentTree[chr];
            vector<Interval*> overlapping;
            tree.Find(start, overlapping);

            for (size_t i = 0; i < overlapping.size(); ++i)
            {
                ExtraInfo& extra = *((ExtraInfo*)(overlapping[i]->mpData));
                extra.mCount++;
                if (strand == ((*extra.mGenes.begin())->mStrand))
                    mSenseReadCount++;
                else
                    mAntiSenseReadCount++;
                //extra.mMapPos.push_back(start - overlapping[i]->mStart);
                //extra.mMapFlag.push_back(flag);
            }
        }
        
    } // OnAShortRead

    void FinishLoadingReads()
    {
        // All the information has been stored in mAllJunctions and mAllSegments
        mSegmentTree.clear();
        mJunctionTree.clear();

        if (mVerbosLevel > 0)
        {
            cerr << "-VL1 " << mSenseReadCount << " sense reads have been counted." << endl;
            cerr << "-VL1 " << mAntiSenseReadCount << " anti-sense reads have been counted." << endl;
        }
    } // FinishLoadingReads

    void BeginLoadingRef()
    {
        mLastChr = "";
        mSegCnt = 0;
        mJuncCnt = 0;
        mLastBedIdx = 0;
        mTailList.clear();
        mAllBedList.clear();

        mpOutput = new ofstream(mOutputFile.data(), ios::out);
        if (!mpOutput->is_open())
        {
            cerr << "File " << mOutputFile << " can not be opened" << endl;
            exit(1);
        }
    }

    void EndLoadingRef()
    {
        mpOutput->close();
        delete mpOutput;
    }

    void OnARefLine(const string& chr, const string& ref_line)
    {
        if (chr != mLastChr)
        {
            //cout << "mAllBedList.size() = " << mAllBedList.size() 
            //     << "\tmLastBedIdx = " << mLastBedIdx << endl;

            assert(mTailList.size() == 0);
            assert((int)mAllBedList.size() == mLastBedIdx);

            mLastChr = chr;

            mAllBedList.resize(0);
            mAllBedList.reserve(mAllSegments[chr].size() + 2 * mAllJunctions[chr].size());
            for_each_ele_in_group(iter, list<Interval>, mAllSegments[chr])
            {
                ExtraInfo& extra = *((ExtraInfo*)(iter->mpData));
                extra.mSecondPart = -1;            // An indicator for segment
                mAllBedList.push_back(*iter);
            }
            for_each_ele_in_group(iter, list<Interval>, mAllJunctions[chr])
            {
                ExtraInfo& extra = *((ExtraInfo*)(iter->mpData));
                mAllBedList.push_back(*iter);
                mAllBedList.push_back(Interval(extra.mSecondPart, extra.mSecondPart + mFlankLen, iter->mpData));
                extra.mSecondPart = -2;            // An indicator for junction waiting for the first part being extracted
            }
            sort(mAllBedList.begin(), mAllBedList.end());

            // Get the length of the longest bed
            mRefSeqBufLen = 0;
            for (size_t i = 0; i < mAllBedList.size(); ++i)
                if (mRefSeqBufLen < mAllBedList[i].mEnd - mAllBedList[i].mStart)
                    mRefSeqBufLen = mAllBedList[i].mEnd - mAllBedList[i].mStart;
            mRefSeqBufLen += ref_line.length() * 10;   // Extend the size a little bit
            mRefSeqBuf.resize(mRefSeqBufLen);
            mLongestLine = ref_line.length();
            mLastBedIdx = 0;
            mLastPos = 0;
        }

        if (mLongestLine < (int)ref_line.length())
            cerr << "WARNING: Your refseq file looks strange, the length of each line is not equal. The result may not be correct." << endl;

        // Save refseq to buffer
        for (size_t i = 0; i < ref_line.length(); ++i)
            mRefSeqBuf[(mLastPos + i) % mRefSeqBufLen] = ref_line[i];

        // Remember intervals with start position in current line
        for (; mLastBedIdx < (int)mAllBedList.size(); ++mLastBedIdx)
        {
            Interval& head_inter = mAllBedList[mLastBedIdx];
            if (mLastPos <= head_inter.mStart && head_inter.mStart < (int)(mLastPos + ref_line.length()))
            {
                // Sort mTailList by the end position
                list<Interval>::iterator iter = mTailList.begin();
                while (iter != mTailList.end() && head_inter.mEnd > iter->mEnd) ++iter;
                mTailList.insert(iter, head_inter);
            }
            else
                break;
        }

        // Get refseq for intervals with end position before the end position of current line
        list<Interval>::iterator iter = mTailList.begin();
        while (iter != mTailList.end() && iter->mEnd <= (int)(mLastPos + ref_line.length()))
        {
            // A ref seq should be extracted.
            string partial_ref_seq;
            partial_ref_seq.resize(iter->mEnd - iter->mStart);
            for (int i = iter->mStart; i < iter->mEnd; ++i)
                partial_ref_seq[i - iter->mStart] = mRefSeqBuf[i % mRefSeqBufLen];

            ExtraInfo& extra = *((ExtraInfo*)(iter->mpData));
            if (-1 == extra.mSecondPart)   // A segment
            {
                OutputASegment(++mSegCnt, *iter, partial_ref_seq);
                partial_ref_seq.clear();
            }
            else if (-2 == extra.mSecondPart)   // A junction waiting for the first part being extracted
            {
                extra.mFirstPartSeq = partial_ref_seq;
                extra.mSecondPart = iter->mEnd;   // Remember the end position of the first part
            }
            else
            {
                OutputAJunction(++mJuncCnt, *iter, partial_ref_seq);
            }
            ++iter;
        }
        mTailList.erase(mTailList.begin(), iter);
        mLastPos += ref_line.length();
    }

    void OutputASegment(int seg_cnt, Interval& inter, string& partial_ref_seq)
    {
        ExtraInfo& extra = *((ExtraInfo*)(inter.mpData)) ;

        list<Gene*>& genes = extra.mGenes;
        list<int>& map_pos = extra.mMapPos;
        list<int>& map_flag = extra.mMapFlag;
        (*mpOutput) << (*genes.begin())->mChr << "\t"  
             << inter.mStart << "\t" 
             << inter.mEnd  << "\t" 
             << "segment" << seg_cnt << "\t"
             << partial_ref_seq << "\t"
             << (*genes.begin())->mStrand << "\t";
        for_each_ele_in_group(gene_iter, list<Gene*>, genes)
            (*mpOutput) << (*gene_iter)->mGeneSymbol<< ",";
        (*mpOutput) << "\t";
        for_each_ele_in_group(gene_iter, list<int>, map_pos)
            (*mpOutput) << (*gene_iter) << ",";
        (*mpOutput) << "\t";
        for_each_ele_in_group(gene_iter, list<int>, map_flag)
            (*mpOutput) << (*gene_iter) << ",";
        (*mpOutput) << endl;
    }

    void OutputAJunction(int junc_cnt, Interval& inter, string& partial_ref_seq)
    {
        ExtraInfo& extra = *((ExtraInfo*)(inter.mpData)) ;
        int first_part_end = extra.mSecondPart;
        int second_part_start = inter.mStart;

        list<Gene*>& genes = extra.mGenes;
        list<int>& map_pos = extra.mMapPos;
        list<int>& map_flag = extra.mMapFlag;
        (*mpOutput) << (*genes.begin())->mChr << "\t"  
             << first_part_end << "\t" 
             << second_part_start << "\t" 
             << "junction" << junc_cnt << "\t"
             << extra.mFirstPartSeq << partial_ref_seq << "\t"
             << (*genes.begin())->mStrand << "\t";
        for_each_ele_in_group(gene_iter, list<Gene*>, genes)
            (*mpOutput) << (*gene_iter)->mGeneSymbol<< ",";
        (*mpOutput) << "\t";
        for_each_ele_in_group(gene_iter, list<int>, map_pos)
            (*mpOutput) << (*gene_iter) << ",";
        (*mpOutput) << "\t";
        for_each_ele_in_group(gene_iter, list<int>, map_flag)
            (*mpOutput) << (*gene_iter) << ",";
        (*mpOutput) << endl;
    }

    // After GPF and SAM are loaded, call this function to print out the number of
    // reads mapped to each gene
    void PrintGeneReadCount()
    {
        typedef map<string, int> str2int_t;
        str2int_t gene_cnt;
        str2int_t gene_length;

        // Note that no need to collect information from mAllJunctions if mbCountJunc is not set
        assert(!mbCountJunc);
        for_each_ele_in_group(iter, map_str2list_inter_t, mAllSegments)
        {
            for_each_ele_in_group(sub_iter, list<Interval>, iter->second)
            {
                ExtraInfo& extra = *((ExtraInfo*)(sub_iter->mpData));
                int cnt = extra.mCount;
                list<Gene*>& genes = extra.mGenes;
                for_each_ele_in_group(gene_iter, list<Gene*>, genes)
                {
                    string& gene_sym = (*gene_iter)->mGeneSymbol;
                    if (gene_cnt.find(gene_sym) == gene_cnt.end())
                    {
                        gene_cnt[gene_sym] = cnt;
                        gene_length[gene_sym] = sub_iter->Length();
                    }
                    else
                    {
                        gene_cnt[gene_sym] += cnt;
                        gene_length[gene_sym] += sub_iter->Length();
                    }
                }
            }
        }

        ofstream output(mOutputFile.data(), ios::out);
        if (!output.is_open())
        {
            cerr << "File " << mOutputFile << " can not be opened" << endl;
            exit(1);
        }

        for_each_ele_in_group(iter, str2int_t, gene_cnt)
        output << iter->first << "\t" << iter->second << "\t" << gene_length[iter->first] << endl;

        output.close();
    }
};

#endif
