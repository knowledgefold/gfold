/*
 * =====================================================================================
 *
 *       Filename:  GFOLD.hpp
 *
 *    Description:  The class for estimate differentiall expressed genes
 *
 *        Version:  1.0
 *        Created:  04/19/11 13:41:14
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Jianxing Feng (), jianxing.tongji@gmail.com
 *        Company:  Tongji Univ.
 *
 * =====================================================================================
 */

#ifndef GFOLD_H
#define GFOLD_H

#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip> 
#include <cstdlib>
#include <assert.h>
#include <limits>
#include <set>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include "DataProcessor.hpp"

using namespace std;

typedef map<string, double> str2double_t;


// =====================================================================================
//        Class:  GFOLD
//  Description:  The class for find differentially expressed genes
// =====================================================================================
class GFOLD
{
private:
    size_t mFirstGroupCnt;
    string mNormalizationMethod;

    string2vec_str_t mGeneSegs;
    vector<string> mAllGeneIDs; 
    vector<string> mAllGeneNames; 

    double mSignificantCutoff;
    int mBurnInCount;
    int mSampledCount;
    int mRandomPairCnt;
    bool mbAccurate;

public:
    int mVerbosLevel;
    string mVersion;

public:
    GFOLD(string version, int verbos_level = 0, string normal_method = "Count", 
          int burn_in_count = 100, int sampled_count = 1000, double significant_cutoff = 0.05, 
	  int random_pair_cnt = 20, bool b_accurate = false)
    {
        mNormalizationMethod = normal_method;
        mVerbosLevel = verbos_level;
        mSignificantCutoff = significant_cutoff;
        mBurnInCount = burn_in_count;
        mSampledCount = sampled_count;
        mVersion = version;
        mRandomPairCnt = random_pair_cnt;
	    mbAccurate = b_accurate;
    }

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateAll
    // Description:  On the case of no replications, use posterior Poisson to estimate the
    //               distribution of gene expression. On the case of with replicates, use
    //               lognormal to estimate the distribution of gene expression
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void CalculateAll(const vector<string>& first_group_samples, const vector<string>& second_group_samples, 
                                           string sample_suffix, string genedescfile, string output_file, string output_file_ext)
    {
        mFirstGroupCnt = first_group_samples.size();

        vector<int> gene_length;
        vector<vector<int> > first_group_gene_read_counts;
        vector<vector<int> > second_group_gene_read_counts;

        first_group_gene_read_counts.resize(first_group_samples.size());
        for (size_t i = 0; i < first_group_samples.size(); ++i)
        { 
            string filename = first_group_samples[i] + sample_suffix;
            vector<string> gene_ids;
            vector<string> gene_names;

            vector<vector<int> > dummy;
            first_group_gene_read_counts[i].resize(0);
            loadGeneReadCounts(filename, gene_ids, gene_names, first_group_gene_read_counts[i], gene_length, mVerbosLevel);

            if (mAllGeneIDs.size() > 0 && mAllGeneIDs != gene_ids)
            {
                cerr << "ERROR: The read count file " << filename << " is not in the right format. Please refer to the documentation." << endl;
                exit(0);
            }
            mAllGeneIDs = gene_ids;
            mAllGeneNames = gene_names;
        }

        second_group_gene_read_counts.resize(second_group_samples.size());
        for (size_t i = 0; i < second_group_samples.size(); ++i)
        { 
            string filename = second_group_samples[i] + sample_suffix;
            vector<string> gene_ids;
            vector<string> gene_names;

            vector<vector<int> > dummy;
            second_group_gene_read_counts[i].resize(0);
            loadGeneReadCounts(filename, gene_ids, gene_names, second_group_gene_read_counts[i], gene_length, mVerbosLevel);

            if (mAllGeneIDs != gene_ids)
            {
                cerr << "ERROR: The read count file " << filename << " is not in the right format. Please refer to the documentation." << endl;
                exit(0);
            }
        }

        // Put all prefixes together
        vector<string> all_samples = first_group_samples;
        all_samples.insert(all_samples.end(), second_group_samples.begin(), second_group_samples.end());

        vector<vector<int> > all_gene_read_counts;
        all_gene_read_counts = first_group_gene_read_counts;
        all_gene_read_counts.insert(all_gene_read_counts.end(), second_group_gene_read_counts.begin(), second_group_gene_read_counts.end());

        vector<double> normalize_constants; 
        vector<int> total_cnt; 
        NormalizeConstant(all_gene_read_counts, all_samples, normalize_constants, total_cnt);

        vector<int> first_total_cnts; 
        vector<int> second_total_cnts; 
        vector<double> normalize_constants_first; 
        vector<double> normalize_constants_second; 
        normalize_constants_first.resize(mFirstGroupCnt);
        normalize_constants_second.resize(normalize_constants.size() - mFirstGroupCnt);
        first_total_cnts.resize(mFirstGroupCnt);
        second_total_cnts = first_total_cnts;
        for (size_t i = 0; i < normalize_constants.size(); ++i)
            if (i < mFirstGroupCnt)
            {
                normalize_constants_first[i] = normalize_constants[i];
                first_total_cnts[i] = total_cnt[i];
            }
            else
            {
                normalize_constants_second[i - mFirstGroupCnt] = normalize_constants[i];
                second_total_cnts[i - mFirstGroupCnt] = total_cnt[i];
            }

        vector<double> first_rpkm;
        vector<double> second_rpkm;
        if (gene_length.size() > 0)
            CalculateRPKM(gene_length, first_group_gene_read_counts, second_group_gene_read_counts, 
                          first_total_cnts, second_total_cnts, first_rpkm, second_rpkm);

        vector<double> gfold_value;
        vector<double> logfdc;
        vector<double> fdr;

        fdr.assign(first_group_gene_read_counts[0].size(), 1);
        if (mbAccurate && first_group_gene_read_counts.size() == 1 && second_group_gene_read_counts.size() == 1)
            CalculateAccurateGFOLD(first_group_gene_read_counts[0], second_group_gene_read_counts[0], 
		 		   normalize_constants_first[0], normalize_constants_second[0], gfold_value, logfdc);
        else
            CalculateMultipleReplicates(first_group_gene_read_counts, second_group_gene_read_counts, 
					normalize_constants_first, normalize_constants_second, gfold_value, logfdc, fdr);

        ofstream output(output_file.data(), ios::out);
        if (!output.is_open())
        {
            cerr << "File " << output_file << " cannot be opened" << endl;
            exit(1);
        }

        time_t rawtime;
        struct tm * timeinfo;
        time ( &rawtime );
        timeinfo = localtime ( &rawtime );

        output << "# This file is generated by gfold " << mVersion << " on " << asctime(timeinfo);
        output << "# Normalization constants :" << endl;
        for (size_t i = 0; i < normalize_constants.size(); ++i)
            output << "#    " << all_samples[i] << "\t" 
                    << total_cnt[i] << "\t" << normalize_constants[i] << endl;
        output << "# The GFOLD value could be considered as a reliable log2 fold change." << endl;
        output << "# It is positive/negative if the gene is up/down regulated." << endl;
        output << "# A gene with zero GFOLD value should never be considered as " << endl;
        output << "# differentially expressed. For a comprehensive description of " << endl;
        output << "# GFOLD, please refer to the manual." << endl;
        output << "#GeneSymbol\tGeneName\tGFOLD("<< mSignificantCutoff << ")\tE-FDR\tlog2fdc";
	    if (gene_length.size() > 0)
            output << "\t1stRPKM\t2ndRPKM";
        output << endl;
        for (size_t i = 0; i < mAllGeneIDs.size(); ++i)
        {
            output << mAllGeneIDs[i] << "\t";
            output << mAllGeneNames[i] << "\t";
            output << gfold_value[i] << "\t";
            output << fdr[i] << "\t";
            output << logfdc[i];
            if (gene_length.size() > 0)
                output << "\t" << first_rpkm[i] << "\t" << second_rpkm[i];
            output << endl;
        }
        output.close();

        // The followling output is for outputing normalized read count and associating gene name to gene description
        map<string, string> id2desc;
        if (genedescfile != "")
        {
            ifstream genedesc(genedescfile.data(), ios::in);
            if (!genedesc.is_open())
            {
                cerr << "File " << genedescfile << " cannot be opened" << endl;
                exit(1);
            }
            else
            {
                string line;
                while ( getline (genedesc, line) )
                {
                    vector<string> fields;
                    split(line, '\t', fields);
                    id2desc[fields[1]] = fields[0];
                }
                genedesc.close();
            }
        }

        ofstream outputext(output_file_ext.data(), ios::out);
        if (!outputext.is_open())
        {
            cerr << "File " << output_file_ext << " cannot be opened" << endl;
            exit(1);
        }

        outputext << "# This file is generated by gfold " << mVersion << " on " << asctime(timeinfo);
        outputext << "# Normalization constants :" << endl;
        for (size_t i = 0; i < normalize_constants.size(); ++i)
            outputext << "#    " << all_samples[i] << "\t" 
                    << total_cnt[i] << "\t" << normalize_constants[i] << endl;

        outputext << "# GeneID";
        for (size_t i = 0; i < normalize_constants.size(); ++i)
            outputext << "\t" << all_samples[i];
        outputext << "\tDescription" << endl;

        for (size_t i = 0; i < mAllGeneIDs.size(); ++i)
        {
            outputext << mAllGeneIDs[i];
            for (size_t j = 0; j < all_gene_read_counts.size(); ++j)
                outputext << "\t" << all_gene_read_counts[j][i] / normalize_constants[j];
            if (id2desc.find(mAllGeneIDs[i]) == id2desc.end())
                outputext << "\tNA";
            else
                outputext << "\t" << id2desc[mAllGeneIDs[i]];
            outputext << endl;
        }
        outputext.close();

    } // CalculateAll

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateSingleReplicate
    // Description:  For the case of no replicates
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void CalculateAccurateGFOLD(const vector<int>& first_group_gene_read_counts, 
				const vector<int>& second_group_gene_read_counts,
  			        double first_normalize_constant, double second_normalize_constant,
				vector<double>& accurate_gfold_value, vector<double>& logfdc)
    {
        if (mVerbosLevel > 0)
            cerr << "-VL1 Calculate accurate GFOLD value ..." << endl;

        accurate_gfold_value.assign(first_group_gene_read_counts.size(), 0);
	AccurateGFOLD(first_group_gene_read_counts, second_group_gene_read_counts, first_normalize_constant, second_normalize_constant, accurate_gfold_value);

	logfdc = accurate_gfold_value;
	for (unsigned i = 0; i < logfdc.size(); ++i)
	    logfdc[i] = log2((second_group_gene_read_counts[i]+1)/second_normalize_constant) - log2((first_group_gene_read_counts[i]+1)/first_normalize_constant);
    } // CalculateSingleReplicate


    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateSingleReplicate
    // Description:  For the case of no replicates
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void CalculateRPKM(const vector<int>& gene_length,
	  	       const vector<vector<int> >& first_group_gene_read_counts, 
		       const vector<vector<int> >& second_group_gene_read_counts,
		       const vector<int>& first_group_total_counts,
		       const vector<int>& second_group_total_counts,
		       vector<double>& first_rpkm, vector<double>& second_rpkm)
    {
        if (mVerbosLevel > 0)
            cerr << "-VL1 Calculate RPKM ..." << endl;

        double first_total_sum = 1;
        double second_total_sum = 1;
        for (unsigned i = 0; i < first_group_total_counts.size(); ++i)
            first_total_sum += first_group_total_counts[i];
        for (unsigned i = 0; i < second_group_total_counts.size(); ++i)
            second_total_sum += second_group_total_counts[i];
        first_total_sum /= 1000000;
        second_total_sum /= 1000000;

        first_rpkm.assign(gene_length.size(), 0);
        second_rpkm = first_rpkm;
        for (unsigned i = 0; i < gene_length.size(); ++i)
        {
            double first_sum = 0;
            for (unsigned j = 0; j < first_group_gene_read_counts.size(); ++j)
                first_sum += first_group_gene_read_counts[j][i];
            first_rpkm[i] = first_sum * 1000 / gene_length[i] / first_total_sum;
            double second_sum = 0;
            for (unsigned j = 0; j < second_group_gene_read_counts.size(); ++j)
                second_sum += second_group_gene_read_counts[j][i];
            second_rpkm[i] = second_sum * 1000 / gene_length[i] / second_total_sum;
        }
    } // CalculateRPKM

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateSingleReplicate
    // Description:  For the case of no replicates
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void CalculateMultipleReplicates(const vector<vector<int> >& first_group_gene_read_counts, 
	  			     const vector<vector<int> >& second_group_gene_read_counts,
 				     const vector<double>& first_normalize_constants, 
				     const vector<double>& second_normalize_constants,
				     vector<double>& gfold_value, vector<double>& logfdc, vector<double>& fdr)
    {

        vector<vector<double> > sampled_log_expression_first; 
        vector<vector<double> > sampled_log_expression_second;

        if (mVerbosLevel > 0)
            cerr << "-VL1 Sampling posterior distribution of log2 fold change ..." << endl;

        double logvar_first = -1;
        double logvar_second = -1;
        if (first_group_gene_read_counts.size() == 1)
            EstimateMCMCPostPoisson(first_group_gene_read_counts[0], first_normalize_constants[0], sampled_log_expression_first);
        else
            EstimateMCMCLogNormal(first_group_gene_read_counts, first_normalize_constants, sampled_log_expression_first, logvar_first);

        if (second_group_gene_read_counts.size() == 1)
            EstimateMCMCPostPoisson(second_group_gene_read_counts[0], second_normalize_constants[0], sampled_log_expression_second);
        else
            EstimateMCMCLogNormal(second_group_gene_read_counts, second_normalize_constants, sampled_log_expression_second, logvar_second);

        if (logvar_first < 0 && logvar_second > 0)
            logvar_first = logvar_second;
        if (logvar_second < 0 && logvar_first > 0)
            logvar_second = logvar_first;

        vector<double> logfdc_low;
        vector<double> logfdc_high;
        vector<double> logfdc_skewness;
        vector<double> logfdc_sd;

        if (mVerbosLevel > 0)
            cerr << "-VL1 Calculating GFOLD value ..." << endl;

        CalculateGFOLD(sampled_log_expression_first, sampled_log_expression_second, 
                       gfold_value, logfdc, logfdc_low, logfdc_high, logfdc_skewness, logfdc_sd);

	fdr.resize(gfold_value.size());
        if (first_group_gene_read_counts.size() > 1 || second_group_gene_read_counts.size() > 1)
        {
            if (mVerbosLevel > 0)
                cerr << "-VL1 Calculating FDR ..." << endl;

            CalculateFDR(first_group_gene_read_counts, second_group_gene_read_counts, 
                         first_normalize_constants, second_normalize_constants,
                         gfold_value, fdr);
        }
    } // CalculateMultipleReplicates


private:

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  NormalizeConstant 
    // Description:  Given two GNBOptimizers which contain all the needed information, 
    //               calculate p-value of differentially expressed for each gene under 
    //               different conditions
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void NormalizeConstant(const vector<vector<int> >& gene_grouped_counts, const vector<string>& samples,
                           vector<double>& normalize_constants, vector<int>& total_cnt)
    {
        if (mVerbosLevel > 0)
            cerr << "-VL1 Calculating normalization constant ..." << endl;

        normalize_constants.resize(gene_grouped_counts.size());
        total_cnt.resize(normalize_constants.size());

        // Always normalize by count
        for (size_t i = 0; i < gene_grouped_counts.size(); ++i)
        {
            int mean = 0;
            for (size_t j = 0; j < gene_grouped_counts[i].size(); ++j)
                mean += gene_grouped_counts[i][j];
            normalize_constants[i] = mean;
            total_cnt[i] = mean;
        }

        if (mNormalizationMethod == "NO")
        {
            normalize_constants.assign(normalize_constants.size(), 1.0);
        }
        else if (mNormalizationMethod == "DESeq")
        {
            bool b_succ = true;
            vector<double> curr_normalize_constants = normalize_constants;
            vector<double> per_gene;
            per_gene.resize(gene_grouped_counts[0].size());
            for (size_t i = 0; i < gene_grouped_counts.size(); ++i)
            {
                for (size_t j = 0; j < gene_grouped_counts[i].size(); ++j)
                {
                    double mean = 0;
                    for (size_t k = 0; k < gene_grouped_counts.size(); ++k)
                        mean += log(gene_grouped_counts[k][j] + 1); // Add a pseudo count
                    mean *= (double)1 / gene_grouped_counts.size();
                    per_gene[j] = gene_grouped_counts[i][j] / exp(mean);
                }
                curr_normalize_constants[i] = UtilityTemp<double>::median(per_gene);

                if (0 == curr_normalize_constants[i])
                {
                    b_succ = false;
                    break; 
                }
            }

            if (b_succ)
                normalize_constants = curr_normalize_constants;
            else
                cerr << "-VL1 WARNING: The median of gene read count is zero. The normalization method proposed"
                     << " by DESeq cannot be used. Use total read count to do normalization." << endl;
        }
        else if (mNormalizationMethod == "TMM")
        {
            const vector<int>& ref_gene_cnt = gene_grouped_counts[0];

            // Treat the first sample as the reference
            int total_ref_cnt = 0;
            for (size_t j = 0; j < ref_gene_cnt.size(); ++j)
                total_ref_cnt += ref_gene_cnt[j];

            for (size_t i = 0; i < gene_grouped_counts.size(); ++i)
            {
                const vector<int>& curr_genes = gene_grouped_counts[i];

                int total_cnt = 0;
                for (size_t j = 0; j < curr_genes.size(); ++j)
                    total_cnt += curr_genes[j];

                vector<double> all_counts;
                vector<double> Mvalue, Avalue;
                Mvalue.reserve(curr_genes.size());
                Avalue.reserve(curr_genes.size());
                for (size_t j = 0; j < curr_genes.size(); ++j)
                {
                    int ref_cnt = ref_gene_cnt[j];
                    int curr_cnt = curr_genes[j];
                    if (0 == ref_cnt || 0 == curr_cnt)
                        continue;
                    double M = (log2(curr_cnt) - log2(total_cnt)) / (log2(ref_cnt) - log2(total_ref_cnt));
                    double A = 0.5 * log2(curr_cnt/total_cnt * ref_cnt/total_ref_cnt);
                    Mvalue.push_back(M);
                    Avalue.push_back(A);
                }
                sort(Mvalue.begin(), Mvalue.end());
                sort(Avalue.begin(), Avalue.end());
                double bot_M_cutoff = Mvalue[size_t(Mvalue.size() * 0.3)];
                double top_M_cutoff = Mvalue[size_t(Mvalue.size() * 0.7)];
                double bot_A_cutoff = Avalue[size_t(Avalue.size() * 0.05)];
                double top_A_cutoff = Avalue[size_t(Avalue.size() * 0.95)];

                double nor = 0;
                double denor = 0;
                for (size_t j = 0; j < curr_genes.size(); ++j)
                {
                    int ref_cnt = ref_gene_cnt[j];
                    int curr_cnt = curr_genes[j];
                    if (0 == ref_cnt || 0 == curr_cnt)
                        continue;
                    double W = (total_cnt - curr_cnt) / (total_cnt * curr_cnt) + (total_ref_cnt - ref_cnt) / (total_ref_cnt * ref_cnt);
                    double M = (log2(curr_cnt) - log2(total_cnt)) / (log2(ref_cnt) - log2(total_ref_cnt));
                    if (M < bot_M_cutoff || M > top_M_cutoff)
                        continue;
                    double A = 0.5 * log2(curr_cnt/total_cnt * ref_cnt/total_ref_cnt);
                    if (A < bot_A_cutoff || A > top_A_cutoff)
                        continue;
                    nor += W * M;
                    denor += W;
                }
                normalize_constants[i] = nor / denor;
            }
        }
        else{
            // The constants are specified explicitly.
            vector<string> normconst;
            split(mNormalizationMethod, ',', normconst);
            for (size_t i = 0; i < normalize_constants.size(); ++i)
                normalize_constants[i] = atof(normconst[i].data());
        }

        double smallest = normalize_constants[0];
        for (size_t i = 1; i < normalize_constants.size(); ++i)
            if (smallest > normalize_constants[i])
                smallest = normalize_constants[i];
        for (size_t i = 0; i < normalize_constants.size(); ++i)
            normalize_constants[i] /= smallest;
        
        if (mVerbosLevel > 0)
        {
            cerr << "-VL1 Normalization constant is: " << endl;
            for (size_t i = 0; i < gene_grouped_counts.size(); ++i)
            {
                double count= 0;
                for (size_t j = 0; j < gene_grouped_counts[i].size(); ++j)
                    count += gene_grouped_counts[i][j];
                cerr << "     " << samples[i] << "    " << count << "    " << normalize_constants[i] << endl;
            }
        }
    } // NormalizeConstant

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  MultiTestCorrection
    // Description:  Do Benjamini, Benjamini and Hochberg multiple testing correction
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void MultiTestCorrection(const vector<double>& ori_pvalue, vector<double>& bh_corrected_pvalue, vector<double>& ben_corrected_pvalue)
    {
        // TODO: This function seems having some bugs
        vector<int> sortedIndex;
        bh_corrected_pvalue = ori_pvalue;
        ben_corrected_pvalue = ori_pvalue;
        UtilityTempComp<double>::Sort(bh_corrected_pvalue, sortedIndex);
        //The assertion fails on two samples exactly the same
        //assert(bh_corrected_pvalue[0] >= bh_corrected_pvalue[bh_corrected_pvalue.size()-1]);
        double total_cnt = (double)sortedIndex.size();
        // Make sure that ori_pvalue[sortedIndex[0]] is the largest pvalue
        double cummin = 1;
        for (size_t i = 0; i < sortedIndex.size(); ++i)
        {
            int idx = sortedIndex[i];
            ben_corrected_pvalue[i] *= total_cnt;
            bh_corrected_pvalue[idx] = ori_pvalue[idx] * total_cnt / (total_cnt - i);
            if (cummin >= bh_corrected_pvalue[idx])
                cummin = bh_corrected_pvalue[idx];
            else
                bh_corrected_pvalue[idx] = cummin;
            if (ben_corrected_pvalue[i] > 1)
                ben_corrected_pvalue[i] = 1;
            if (bh_corrected_pvalue[sortedIndex[i]] > 1)
                bh_corrected_pvalue[sortedIndex[i]] = 1;
        }
    } // MultiTestCorrection


    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateGFOLD
    // Description:  Given two set of sampled posterior distribution, calculate the GFOLD value
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void CalculateGFOLD(const vector<vector<double> >& sampled_log_expression_first,
                       const vector<vector<double> >& sampled_log_expression_second,
                       vector<double>& GFOLD,
                       vector<double>& logfdc,
                       vector<double>& logfdc_low,
                       vector<double>& logfdc_high,
                       vector<double>& logfdc_skewness,
                       vector<double>& logfdc_sd)
    {
        size_t sampled_cnt = sampled_log_expression_first[0].size();

        logfdc.resize(mAllGeneIDs.size());
        logfdc_sd = logfdc_skewness = logfdc_low = logfdc_high = GFOLD = logfdc;
        for (size_t i = 0; i < mAllGeneIDs.size(); ++i)
        {
            vector<double> lfdc;
            lfdc.resize(sampled_cnt);

            double values[sampled_cnt];
            for (size_t j = 0; j < sampled_cnt; ++j)
            {
                lfdc[j] = sampled_log_expression_second[i][j] - sampled_log_expression_first[i][j];
                lfdc[j] /= log(2);
                values[j] = lfdc[j];
            }

            logfdc[i] = UtilityTemp<double>::mean(lfdc);
            sort(lfdc.begin(), lfdc.end());
            logfdc_low[i] = lfdc[size_t((lfdc.size()-1) * mSignificantCutoff)];
            logfdc_high[i] = lfdc[size_t((lfdc.size()-1) * (1-mSignificantCutoff))];

            logfdc_skewness[i] = gsl_stats_skew(values, 1, sampled_cnt);
            logfdc_sd[i] = gsl_stats_sd(values, 1, sampled_cnt);

            if (0 >= logfdc_low[i] && 0 <= logfdc_high[i] ) 
                GFOLD[i] = 0;
            else if (logfdc_low[i] > 0)
                GFOLD[i] = logfdc_low[i];
            else
                GFOLD[i] = logfdc_high[i];
        }
    }

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  LogFactorial
    // Description:  Approximate the log factorial
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    double LogFactorial(int n)
    {
        double value;
        if (n < 10)
        {
            value = 0;
            for (int i = 1; i < n; ++i)
               value += log(i); 
        }
        else
            value = n * log(n) - n + (log(n) + log(4*n) + log(1+2*n))/6 + log(3.1415926) / 2;
        return value;
    } //LogFactorial


    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  AccurateProb
    // Description:  Given a cutoff, raw read counts and normalization constant
    //               calculate the probability
    //               P(\lambda1 \geq \alpha*\lambda2)
    //               where \lambda1 is the posterior distribution of lambda under poission distribution
    //               if the observed read count is k1 and the normalization constant is nc1.
    //
    //               Let \lambda1' = \lambda1 * nc1
    //               Let \lambda2' = \lambda2 * nc2
    //                 P(\lambda1 \geq \alpha*\lambda2)
    //               = P(\lambda1'/nc1 \geq \alpha*\lambda2'/nc2)
    //               = P(\lambda1' \geq \alpha*nc1/nc2*\lambda2')
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    double AccurateProb(int k1, int k2, double alpha, double nc1, double nc2)
    {
        double prob = 0;
        alpha *= nc1 / nc2;

        if (k1 < k2)
        {
            for (int i = 0; i <= k1; ++i)
            {
                double curr = i * log(alpha/(alpha+1)) + LogFactorial(k2+i) - LogFactorial(i) - LogFactorial(k2) - (k2+1) * log(alpha+1);
                //cerr << k1 << "\t" << k2 << "\t" << i << "\t" << alpha << "\t" << curr << endl;
                prob += exp(curr); 
            }
        }
        else
        {
            // Use the fact that P(\lambda1 \geq \alpha*\lambda2) = 1 - P(\lambda2 \geq 1/\alpha*\lambda1)
            alpha = 1 / alpha;
            for (int i = 0; i <= k2; ++i)
            {
                double curr = i * log(alpha/(alpha+1)) + LogFactorial(k1+i) - LogFactorial(i) - LogFactorial(k1) - (k1+1) * log(alpha+1);
                prob += exp(curr); 
            }

            prob = 1 - prob;
        }

        //cerr << k1 << "\t" << k2 << "\t" << alpha << "\t" << prob << endl;

        return prob;
    } //AccurateProb



    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  AccurateGFOLD
    // Description:  Calculate the accurate gfold value
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    void AccurateGFOLD(const vector<int>& sample1_raw_read_counts, 
                      const vector<int>& sample2_raw_read_counts, 
                      double nc1, 
                      double nc2, 
                      vector<double>& accurate_gfold_value)
    {
        accurate_gfold_value.assign(sample1_raw_read_counts.size(), 0);

        for (size_t i = 0; i < sample1_raw_read_counts.size(); ++i)
        {
            int cnt1 = sample1_raw_read_counts[i];
            int cnt2 = sample2_raw_read_counts[i];
            double low_cutoff = 0;
            double precision = 20;
            double high_cutoff = 0;
            double alpha;

            // Down regulated
            if (cnt1 / nc1 > cnt2 / nc2)
            {
                if (cnt2 == 0) cnt2 += 1;
                high_cutoff = (cnt1/nc1/(cnt2/nc2));
                for (int j = 0; j < precision; ++j)
                {
                    alpha = (low_cutoff + high_cutoff) / 2;
                    if (1 - mSignificantCutoff < AccurateProb(cnt1, cnt2, alpha, nc1, nc2))
                        low_cutoff = alpha;
                    else
                        high_cutoff = alpha;
                }

                accurate_gfold_value[i] = -log2(alpha);
                if (accurate_gfold_value[i] > 0)
                    accurate_gfold_value[i] = 0;
            }
            else // Up regulated
            {
                if (cnt1 == 0) cnt1 += 1;
                high_cutoff = (cnt2/nc2/(cnt1/nc1));
                for (int j = 0; j < precision; ++j)
                {
                    alpha = (low_cutoff + high_cutoff) / 2;
                    if (1 - mSignificantCutoff < AccurateProb(cnt2, cnt1, alpha, nc2, nc1))
                        low_cutoff = alpha;
                    else
                        high_cutoff = alpha;
                }

                accurate_gfold_value[i] = log2(alpha);
                if (accurate_gfold_value[i] < 0)
                    accurate_gfold_value[i] = 0;
            }
        }
    } //AccurateGFOLD


    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  EstimateMCMCPostPoisson
    // Description:  If there is no replicate, do posterior Poisson sampling
    //  Parameters:  
    //--------------------------------------------------------------------------------------
    int EstimateMCMCPostPoisson(const vector<int>& raw_read_counts, double NC, 
                                vector<vector<double> >& sampled_log_expression)
    {
        // Remove genes with no read counts at all.
        vector<int> read_counts = raw_read_counts;
        for (size_t i = 0; i < read_counts.size(); ++i)
            read_counts[i] += 1;

        vector<double> expression;
        expression.resize(read_counts.size());

        vector<double> proposed_sd_expression = expression;

        for (size_t i = 0; i < read_counts.size(); ++i)
        {
            expression[i] = (read_counts[i] + 1) / NC;
            proposed_sd_expression[i] = sqrt(expression[i]);
        }

        int sampled_count = 1000;

        vector<vector<double> > saved;
        saved.resize(sampled_count);

        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        double learn_step_inc = 0.1;
        double learn_step_dec = 0.1;
        int reject_cnt = 0;
        int iter = 0;
        while (iter++ < sampled_count)
        {
            for (size_t i = 0; i < expression.size(); ++i)
            {
                double old_value = expression[i];
                double proposed_mean = old_value; 
                double proposed_std = proposed_sd_expression[i]; 
                double new_value = -1;
                while (new_value <= 0)
                    new_value = gsl_ran_gaussian(r, 1) * proposed_std + proposed_mean;

                // Weight for IS
                double old_log_prob = pow((old_value - proposed_mean) / proposed_std, 2) / 2;
                double new_log_prob = pow((new_value - proposed_mean) / proposed_std, 2) / 2;
                // Prior is uniform
                // Likelihood
                old_log_prob += read_counts[i] * log(old_value * NC) - old_value * NC;
                new_log_prob += read_counts[i] * log(new_value * NC) - new_value * NC;

                double u = gsl_rng_uniform (r);
                if (u < exp(new_log_prob - old_log_prob))
                {
                    expression[i] = new_value;
                    proposed_sd_expression[i] *= 1 + learn_step_inc;
                }
                else
                {
                    ++reject_cnt;
                    proposed_sd_expression[i] *= 1 - learn_step_dec;
                }
            }

            // Save sampled elements
            int curr_idx = iter % sampled_count;
            saved[curr_idx] = expression;

            if (mVerbosLevel > 2 && iter % 10 == 0)
            {
                cerr << "---VL3 iter = " << iter 
                     << "  NC = " << NC
                     << "  expression[10000] = " << expression[10000] << endl;
            }
        }

        sampled_log_expression.resize(read_counts.size());

        for (size_t i = 0; i < sampled_log_expression.size(); ++i) 
        {
            sampled_log_expression[i].resize(sampled_count);
            for (int j = 0; j < sampled_count; ++j) 
                sampled_log_expression[i][j] = log(saved[j][i]);
        }

        gsl_rng_free (r);

        return sampled_count;
    } //EstimateMCMCPostPoisson

    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  EstimateMCMCLogNormal
    // Description:  Estimate parameters using MCMC. This method use data augmentation
    //               algorithm with true expression of each gene as the augmented parameter (the
    //               mission values).  In each round, the posterior is calculated using Metropolized 
    //               independence sampler (MIS). More information could be found in "Monte Carlo 
    //               Strategies in Scientific Computing" by Jun S. Liu, Page 115, Page 136, Page 312.
    //
    //  Parameters:  [in] raw_read_counts: The read count of each element in exp_or_bias
    //               [in] normalization_constant: The normalization constant for each sample
    //               [out] sampled_log_expression_mean: The sampled expression of each gene in log scale when both sampled are considered
    //
    //      Return:  iterate times
    //--------------------------------------------------------------------------------------
    int EstimateMCMCLogNormal(const vector<vector<int> >& raw_read_counts, 
                              const vector<double>& normalization_constant, 
                              vector<vector<double> >& sampled_log_expression_mean,
                              double& logvar)
    {
        vector<vector<int> > read_counts = raw_read_counts;
        for (size_t i = 0; i < read_counts.size(); ++i)
            for (size_t j = 0; j < read_counts[i].size(); ++j)
                read_counts[i][j] += 1;

        int gene_cnt = read_counts[0].size();

        // The expression of each gene in each sample
        vector<vector<double> > expression;        
        expression.resize(read_counts.size());
        for (size_t i = 0; i < expression.size(); ++i)
            expression[i].resize(gene_cnt);
        vector<vector<double> > proposed_sd_expression = expression;

        for (size_t i = 0; i < read_counts[0].size(); ++i)
        {
            for (size_t sc = 0; sc < read_counts.size(); ++sc)
            {
                expression[sc][i] = (read_counts[sc][i] + 1) / normalization_constant[sc];
                proposed_sd_expression[sc][i] = sqrt(expression[sc][i]);
            }
        }

        // The mean of the log expression of each gene
        vector<double> logmean;  
        logmean.resize(gene_cnt);
        vector<double> proposed_sd_logmean = logmean;
        for (size_t i = 0; i < read_counts[0].size(); ++i)
        {
            for (size_t sc = 0; sc < read_counts.size(); ++sc)
                logmean[i] += log(expression[0][i]);
            logmean[i] /= read_counts.size();
            proposed_sd_logmean[i] = 1;
        }

        // The variance of the log expression of each gene
        double proposed_sd_logvar = 1;

        // The mean of the log expression of all genes
        double proposed_sd_mean_root = 2;

        // The variance of the log expression of all genes
        double proposed_sd_var_root = 2;
        
        logvar = 0.20;
        double mean_root = 4;
        double var_root = 2;
        bool b_draw_logvar = true;

        int burn_in_count = mBurnInCount;
        int sampled_count = mSampledCount;

        vector<vector<double> > sampled_logmean;
        sampled_logmean.resize(sampled_count);

        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);

        double learn_step_inc = 0.1;
        double learn_step_dec = 0.1;
        int reject_cnt = 0;
        int iter = 0;
        while (iter++ < burn_in_count + sampled_count)
        {
            // Draw expression
            for (size_t sc = 0; sc < read_counts.size(); ++sc)
            {
                double norm = normalization_constant[sc];
                for (size_t i = 0; i < expression[sc].size(); ++i)
                {
                    double old_value = expression[sc][i];
                    double proposed_mean = old_value; 
                    double proposed_std = proposed_sd_expression[sc][i]; 
                    double new_value = -1;
                    while (new_value <= 0)
                        new_value = gsl_ran_gaussian(r, 1) * proposed_std + proposed_mean;

                    // expression1[i] ~ dlnorm(logmean[i], revvar[i])
                    // Weight for MIS
                    double old_log_prob = pow((old_value - proposed_mean) / proposed_std, 2) / 2;
                    double new_log_prob = pow((new_value - proposed_mean) / proposed_std, 2) / 2;
                    // Prior
                    old_log_prob += -log(old_value) - pow(log(old_value) - logmean[i], 2) / (2 * logvar);
                    new_log_prob += -log(new_value) - pow(log(new_value) - logmean[i], 2) / (2 * logvar);
                    // Likelihood
                    old_log_prob += read_counts[sc][i] * log(old_value * norm) - old_value * norm;
                    new_log_prob += read_counts[sc][i] * log(new_value * norm) - new_value * norm;

                    double u = gsl_rng_uniform (r);
                    if (u < exp(new_log_prob - old_log_prob))
                    {
                        expression[sc][i] = new_value;
                        proposed_sd_expression[sc][i] *= 1 + learn_step_inc;
                    }
                    else
                    {
                        ++reject_cnt;
                        proposed_sd_expression[sc][i] *= 1 - learn_step_dec;
                    }
                }
            }
  
            // Draw logmean
            for (size_t i = 0; i < logmean.size(); ++i)
            {
                double old_value = logmean[i];
                double proposed_mean = logmean[i];
                double proposed_std = proposed_sd_logmean[i]; 
                double new_value = gsl_ran_gaussian(r, 1) * proposed_std + proposed_mean;

                // expression1[i] ~ dlnorm(logmean[i], revvar[i])
                // Weight for MIS
                double old_log_prob = pow((old_value - proposed_mean) / proposed_std, 2) / 2;
                double new_log_prob = pow((new_value - proposed_mean) / proposed_std, 2) / 2;
                // Prior
                old_log_prob += -pow(old_value - mean_root, 2) / (2 * var_root);
                new_log_prob += -pow(new_value - mean_root, 2) / (2 * var_root);
                // Likelihood
                for (size_t sc = 0; sc < expression.size(); ++sc)
                {
                    double expr = expression[sc][i];
                    old_log_prob += -log(expr) - 0.5*log(logvar) - pow(log(expr) - old_value, 2) / (2 * logvar);
                    new_log_prob += -log(expr) - 0.5*log(logvar) - pow(log(expr) - new_value, 2) / (2 * logvar);
                }

                double u = gsl_rng_uniform (r);
                if (u < exp(new_log_prob - old_log_prob))
                {
                    logmean[i] = new_value;
                    proposed_sd_logmean[i] *= 1 + learn_step_inc;
                }
                else
                {
                    ++reject_cnt;
                    proposed_sd_logmean[i] *= 1 - learn_step_dec;
                }
            }
 
            // Draw logvar
            if (b_draw_logvar)
            {
                double old_value = logvar;
                double proposed_mean = old_value;
                double proposed_std = proposed_sd_logvar; 
                double new_value = -1;
                while (new_value < 0)
                    new_value = gsl_ran_gaussian(r, 1) * proposed_std + proposed_mean;

                // expression1[i] ~ dlnorm(logmean[i], revvar[i])
                // Weight for MIS
                double old_log_prob = pow((old_value - proposed_mean) / proposed_std, 2) / 2;
                double new_log_prob = pow((new_value - proposed_mean) / proposed_std, 2) / 2;
                // Prior is uniform
                // Likelihood
                for (size_t i = 0; i < expression[0].size(); ++i)
                {
                    for (size_t sc = 0; sc < expression.size(); ++sc)
                    {
                        double expr = expression[sc][i];
                        old_log_prob += -0.5*log(old_value) - pow(log(expr) - logmean[i], 2) / (2 * old_value);
                        new_log_prob += -0.5*log(new_value) - pow(log(expr) - logmean[i], 2) / (2 * new_value);
                    }
                }

                double u = gsl_rng_uniform (r);
                if (u < exp(new_log_prob - old_log_prob))
                {
                    logvar = new_value;
                    proposed_sd_logvar *= 1 + learn_step_inc;
                }
                else
                {
                    ++reject_cnt;
                    proposed_sd_logvar *= 1 - learn_step_dec;
                }
            }
  
            // Draw mean_root
            {
                double old_value = mean_root;
                double proposed_mean = old_value;
                double proposed_std = proposed_sd_mean_root; 
                double new_value = gsl_ran_gaussian(r, 1) * proposed_std + proposed_mean;

                // Weight for MIS
                double old_log_prob = pow((old_value - proposed_mean) / proposed_std, 2) / 2;
                double new_log_prob = pow((new_value - proposed_mean) / proposed_std, 2) / 2;
                // Prior is uniform
                // Likelihood
                for (size_t i = 0; i < logmean.size(); ++i)
                {
                    old_log_prob += - pow(logmean[i] - old_value, 2) / (2 * var_root);
                    new_log_prob += - pow(logmean[i] - new_value, 2) / (2 * var_root);
                }

                double u = gsl_rng_uniform (r);
                if (u < exp(new_log_prob - old_log_prob))
                {
                    mean_root = new_value;
                    proposed_sd_mean_root *= 1 + learn_step_inc;
                }
                else
                {
                    ++reject_cnt;
                    proposed_sd_mean_root *= 1 - learn_step_dec;
                }
            }

            // Draw var_root
            {
                double old_value = var_root;
                double b = proposed_sd_var_root / var_root;
                double a = var_root / b;
                double new_value = gsl_ran_gamma(r, a, b); 

                // Weight for MIS
                double old_log_prob = -(a-1) * log(old_value) + old_value / b;
                double new_log_prob = -(a-1) * log(new_value) + new_value / b;
                // Prior is uniform
                // Likelihood
                for (size_t i = 0; i < logmean.size(); ++i)
                {
                    old_log_prob += - pow(logmean[i] - mean_root, 2) / (2 * old_value) - sqrt(old_value);
                    new_log_prob += - pow(logmean[i] - mean_root, 2) / (2 * new_value) - sqrt(new_value);
                }

                double u = gsl_rng_uniform (r);
                if (u < exp(new_log_prob - old_log_prob))
                {
                    var_root = new_value;
                    proposed_sd_var_root *= 1 + learn_step_inc;
                }
                else
                {
                    ++reject_cnt;
                    proposed_sd_var_root *= 1 - learn_step_dec;
                }
            }


            // Save sampled elements
            int curr_idx = iter % sampled_count;
            sampled_logmean[curr_idx] = logmean;

            if (mVerbosLevel > 2 && iter % 10 == 0)
            {
                int idx = 15524;
                cerr << "---VL3 iter = " << iter 
                     << "  mean_root  = " << mean_root
                     << "  var_root  = " << var_root
                     << "  idx = " << idx 
                     << "  logmean = " << logmean[idx]
                     << "  logvar = " << logvar
                     << "  expression =";
                for (size_t i = 0; i < expression.size(); ++i)
                    cerr << " " << expression[i][idx];
                cerr<< endl;
            }
        }

        sampled_log_expression_mean.resize(read_counts[0].size());
        for (size_t i = 0; i < sampled_log_expression_mean.size(); ++i) 
        {
            sampled_log_expression_mean[i].resize(sampled_count);
            for (int j = 0; j < sampled_count; ++j) 
                sampled_log_expression_mean[i][j] = sampled_logmean[j][i];
        }

        gsl_rng_free (r);

        return burn_in_count + sampled_count;
    } //EstimateMCMCLogNormal


    //--------------------------------------------------------------------------------------
    //       Class:  GFOLD
    //      Method:  CalculateFDR
    // Description:  
    //  Parameters:  [in] raw_read_counts: The read count of each element in exp_or_bias
    //               [in] normalization_constant: The normalization constant for each sample
    //               [out] sampled_log_expression_mean: The sampled expression of each gene in log scale when both sampled are considered
    //
    //      Return:  iterate times
    //--------------------------------------------------------------------------------------
    void CalculateFDR(const vector<vector<int> >& first_group_gene_read_counts, 
                      const vector<vector<int> >& second_group_gene_read_counts, 
                      const vector<double>& normalize_constants_first,
                      const vector<double>& normalize_constants_second,
                      const vector<double>& GFOLD,
                      vector<double>& fdr)
    {
        int first_pair = first_group_gene_read_counts.size() * (first_group_gene_read_counts.size() - 1)/2;
        int second_pair = second_group_gene_read_counts.size() * (second_group_gene_read_counts.size() - 1)/2;
        int total_pair = first_pair + second_pair;

        vector<int> indexes;
        indexes.resize(total_pair);
        for (int i = 1; i <= total_pair; ++i)
            indexes[i-1] = i;
        UtilityTemp<int>::Shuffle(indexes);

        // Take the first 20 pairs
        int rc = mRandomPairCnt;
        if (rc > (int)indexes.size())
            rc = (int)indexes.size();

        vector<double> background_gfold;
        background_gfold.reserve(rc * GFOLD.size());

        for (int k = 0; k < rc; ++k)
        {
            int idx = indexes[k];
            bool b_first = true;
            if (idx > first_pair)
            {
                idx -= first_pair;
                b_first = false;
            }

            int i = int(sqrt(2*idx - 0.75) + 0.5);
            int j = idx - (i * i - i + 1)/2 - 1;


            if (mVerbosLevel > 1)
            {
                cerr << "--VL2 sample pair " << i << "," << j << " in the";
                if (b_first) 
                    cerr << " first";
                else
                    cerr << " second";
                cerr << " group has been selected." << endl;
            }

            vector<vector<double> > sampled_log_expression_first;
            vector<vector<double> > sampled_log_expression_second;
            if (b_first)
            {
                EstimateMCMCPostPoisson(first_group_gene_read_counts[i], normalize_constants_first[i], sampled_log_expression_first);
                EstimateMCMCPostPoisson(first_group_gene_read_counts[j], normalize_constants_first[j], sampled_log_expression_second);
            }
            else
            {
                EstimateMCMCPostPoisson(second_group_gene_read_counts[i], normalize_constants_second[i], sampled_log_expression_first);
                EstimateMCMCPostPoisson(second_group_gene_read_counts[j], normalize_constants_second[j], sampled_log_expression_second);
            }

            vector<double> temp_gfold;
            vector<double> logfdc;
            vector<double> logfdc_low;
            vector<double> logfdc_high;
            vector<double> logfdc_skewness;
            vector<double> logfdc_sd;

            CalculateGFOLD(sampled_log_expression_first, sampled_log_expression_second, 
                           temp_gfold, logfdc, logfdc_low, logfdc_high, logfdc_skewness, logfdc_sd);

            background_gfold.insert(background_gfold.end(), temp_gfold.begin(), temp_gfold.end());
        }

        for (size_t i = 0; i < background_gfold.size(); ++i)
            if (background_gfold[i] < 0)
                background_gfold[i] = -background_gfold[i];

        vector<double> abs_gfold = GFOLD;
        for (size_t i = 0; i < abs_gfold.size(); ++i)
            if (abs_gfold[i] < 0)
                abs_gfold[i] = -abs_gfold[i];

        sort(background_gfold.begin(), background_gfold.end());

        fdr = abs_gfold;
        vector<double> sorted_gfold = abs_gfold;
        sort(sorted_gfold.begin(), sorted_gfold.end());

        for (size_t i = 0; i < GFOLD.size(); ++i)
        {
            double gfold = GFOLD[i];
            if (gfold < 0)
                gfold = -gfold;
            int idx1 = UtilityTemp<double>::BinarySearch(background_gfold, gfold);
            int idx2 = UtilityTemp<double>::BinarySearch(sorted_gfold, gfold);
            double ratio1 = (double)(background_gfold.size() - idx1) / background_gfold.size();
            double ratio2 = (double)(sorted_gfold.size() - idx2) / sorted_gfold.size();

            fdr[i] = ratio1 / ratio2;
        }

        sorted_gfold = abs_gfold;
        vector<int> sortedIndex;
        UtilityTempComp<double>::Sort(sorted_gfold, sortedIndex);

        // For each FDR replace it with the smallest FDR sorted later
        double smallest = 1;
        for (int i = (int)(sortedIndex.size() - 1); i >= 0; --i)
        {
            int gene_idx = sortedIndex[i];
            if (smallest > fdr[gene_idx])
                smallest = fdr[gene_idx];
            fdr[gene_idx] = smallest;
        }
    }
};

#endif 
