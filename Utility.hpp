/*
 * Copyright (C) 2012-2017 Jianxing Feng
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

// =====================================================================================
// 
//       Filename:  Utility.hpp
// 
//    Description:  Some utilities
// 
//        Version:  1.0
//        Created:  04/07/2011 13:00:33
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Jianxing Feng (), jianxing.tongji@gmail.com
//        Company:  Tongji Univ.
// 
// =====================================================================================

#ifndef Utility_H
#define Utility_H

#define for_each_ele_in_group_temp(_Iter, _GroupType, _Group) \
	for (typename _GroupType::iterator _Iter = (_Group).begin(); _Iter != (_Group).end(); _Iter++)

#define for_each_ele_in_group_const_temp(_Iter, _GroupType, _Group) \
	for (typename _GroupType::const_iterator _Iter = (_Group).begin(); _Iter != (_Group).end(); _Iter++)

#define for_each_ele_in_group(_Iter, _GroupType, _Group) \
	for (_GroupType::iterator _Iter = (_Group).begin(); _Iter != (_Group).end(); _Iter++)

#define for_each_ele_in_group_const(_Iter, _GroupType, _Group) \
	for (_GroupType::const_iterator _Iter = (_Group).begin(); _Iter != (_Group).end(); _Iter++)

//#define DEBUG

#ifndef DEBUG
	#define TRACE(_x)
#else
	#define TRACE(_x) cout << _x << endl;
#endif

#ifndef DEBUG
	#define LOCATION(_x)
#else
	#define LOCATION(_x) cout << "LOCATION: " << __FILE__ << ":" << __LINE__ << " " << _x << endl;
#endif

#ifndef DEBUG
	#define LOCATE
#else
	#define LOCATE cout << "LOCATION: " << __FILE__ << ":" << __LINE__ << endl;
#endif

#ifndef DEBUG
	#define ASSERT(_x, _y)
#else
	#define ASSERT(_x, _y) \
		if (!(_x)) cout << _y << endl; \
		assert(_x) 
#endif

#include <string>
#include <vector>
#include <stdlib.h>
#include <cstdlib>
#include <assert.h>
#include <map>
#include <cstdlib>
#include <math.h>
#include <algorithm>
#include <gsl/gsl_statistics_int.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>

using namespace std;

/*
 * Splite a string by a single char
 */
void split(const string& a_string, char sep, vector<string>& fields)
{
    fields.clear();
    if (a_string.length() == 0)
        return;
    size_t last = 0;
    for (size_t i = 0; i < a_string.length(); ++i)
    {
        if (a_string[i] == sep)
        {
            fields.push_back(a_string.substr(last, i - last));
            last = i+1;
        }
    }
    if (last < a_string.length())
        fields.push_back(a_string.substr(last, a_string.length() - last));
    else if (a_string[a_string.length()-1] == sep)
        fields.push_back("");
}

void initPreComputedTable(vector<vector<double> >& pre_computed_table, int taylor_precision, int pre_computed_len)
{    
    pre_computed_table.resize(taylor_precision);
    for (size_t i = 0; i < pre_computed_table.size(); ++i)
        pre_computed_table[i].assign(pre_computed_len, 0);

    // pre_computed_table[0][i] = sum_{j=1}^i log(j)
    for (size_t i = 1; i < pre_computed_table[0].size(); ++i)
        pre_computed_table[0][i] = log(i);

    // pre_computed_table[1][i] = sum_{j=1}^i 1/j
    // pre_computed_table[2][i] = sum_{j=1}^i 1/(2*j^2)
    // pre_computed_table[3][i] = sum_{j=1}^i 1/(3*j^3)
    // ...
    for (size_t i = 1; i < pre_computed_table.size(); ++i)
        for (size_t j = 1; j < pre_computed_table[i].size(); ++j)
            pre_computed_table[i][j] = (double)1/i/pow(j,i);


    for (size_t i = 0; i < pre_computed_table.size(); ++i)
        for (size_t j = 1; j < pre_computed_table[i].size(); ++j)
            pre_computed_table[i][j] += pre_computed_table[i][j-1];
}

void initPreComputedTableFrac(vector<vector<double> >& pre_computed_table_frac, int taylor_precision, int pre_computed_len)
{    
    pre_computed_table_frac.resize(taylor_precision);
    for (size_t i = 0; i < pre_computed_table_frac.size(); ++i)
        pre_computed_table_frac[i].assign(pre_computed_len, 0);

    // pre_computed_table_frac[0][i] = sum_{j=1}^i 1/j
    // pre_computed_table_frac[1][i] = sum_{j=1}^i 1/(j^2)
    // pre_computed_table_frac[2][i] = sum_{j=1}^i 1/(j^3)
    // pre_computed_table_frac[3][i] = sum_{j=1}^i 1/(j^4)
    // ...
    for (size_t i = 0; i < pre_computed_table_frac.size(); ++i)
        for (size_t j = 1; j < pre_computed_table_frac[i].size(); ++j)
        {
            pre_computed_table_frac[i][j] = (double)1/pow(j,i+1);
            pre_computed_table_frac[i][j] += pre_computed_table_frac[i][j-1];
        }
}

/*
 * Calculate \sum_{i=1}^(r-1) log(i+b) using Tayler expansion
 */
double logSum(vector<vector<double> >& pre_computed_table, int r, double b)
{
    if (!(r >= 1 && b >= 0))
    {
        cout << r << "\t" << b << endl;
        assert(r >= 1 && b >= 0);
    }
    if (r < 20)
    {
        double obj = 0;
        for (int j = 1; j < r; ++j)
            obj += log(j+b);
        return obj;
    }

    int int_part = (int)b;    // Round down
    int low = int_part;         
    int high = r - 1 + int_part;   
    double frac = b - int_part;
    if (frac > 0.5)
    {
        frac = b - int_part - 1;
        low = int_part + 1;  
        high = r + int_part;
    }

    if ((size_t)high >= pre_computed_table[0].size())
    {
        // Approximate using Stirling's formula
        double n1 = r + b -1; 
        double n2 = b; 
        double e = 2.718282;
        double pi = 3.1415926;
        double obj = n1 * log(n1 / e) + 0.5 * log(2 * pi * n1) + log(1+1/(12*n1)) - 
                    (n2 * log(n2 / e) + 0.5 * log(2 * pi * n2) + log(1+1/(12*n2))); 
        return obj;

#ifdef DEBUG
        cerr << __func__ << "  Warning: r - 1 + int(b) = " << high 
            << " is too large compared to precomputed table " << pre_computed_table[0].size() - 1 
            << " r = " << r 
            << " b = " << b << endl;
#endif
        //double obj = 0;
        //for (int j = 1; j < r; ++j)
        //    obj += log(j+b);
        ////obj += ;
        ////if (low >= pre_computed_table[0].size())
        ////    obj -= 
        ////else
        ////    obj -= logSum(pre_computed_table, (int)b+1, b - (int)b);
        //return obj;
    }

    double obj = 0;
    double neg = 1;
    obj += (pre_computed_table[0][high] - pre_computed_table[0][low]);
    for (size_t i = 1; i < pre_computed_table.size(); ++i)
    {
        obj += neg * pow(frac, i) * (pre_computed_table[i][high] - pre_computed_table[i][low]);
        neg *= -1;
    }
    return obj;
}

/*
 * Calculate \sum_{i=1}^{r-1} 1 / (j+b) using Tayler expansion
 */
double fracSum(vector<vector<double> >& pre_computed_table, int r, double b)
{
    assert(r >= 1 && b >= 0);
    if (r < 20)
    {
        double obj = 0;
        for (int j = 1; j < r; ++j)
            obj += 1/(j+b);
        return obj;
    }

    int int_part = (int)b;    // Round down
    int low = int_part;         
    int high = r - 1 + int_part;   
    double frac = b - int_part;
    if (frac > 0.5)
    {
        frac = b - int_part - 1;
        low = int_part + 1;  
        high = r + int_part;
    }

    if ((size_t)high >= pre_computed_table[0].size())
    {
        return log(r+b-1) - log(b);

#ifdef DEBUG
        cerr << __func__ << "  Warning: r - 1 + int(b) = " << high 
            << " is too large compared to precomputed table " << pre_computed_table[0].size() - 1 
            << " r = " << r 
            << " b = " << b << endl;
#endif
        double obj = 0;
        for (int j = 1; j < r; ++j)
            obj += 1/(j+b);
        return obj;
    }

    double obj = 0;
    double neg = 1;
    for (size_t i = 0; i < pre_computed_table.size(); ++i)
    {
        obj += neg * pow(frac, i) * (pre_computed_table[i][high] - pre_computed_table[i][low]);
        neg *= -1;
    }
    return obj;
}

// Shuffle all the weight and associated count
void shuffle_weight_count(vector<double>& weight, vector<int>& read_counts)
{
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    for (size_t i = 0; i < weight.size(); ++i)
    {
        double u = gsl_rng_uniform (r);
        int idx = i + (int)(u * (weight.size() - i));

        double temp = weight[i];
        weight[i] = weight[idx];
        weight[idx] = temp;

        int temp_r = read_counts[i];
        read_counts[i] = read_counts[idx];
        read_counts[idx] = temp_r;
    }
    gsl_rng_free (r);
}


template <typename Type>
class UtilityTemp
{
public:
    /*
     * Calculate the sum of a group of values
     */
    inline static 
    Type sum(const vector<Type>& values)
    {
        Type mean = 0;
        for (size_t i = 0; i < values.size(); ++i)
            mean += values[i];
        return mean;
    }

    /*
     * Calculate the mean of a group of values
     */
    inline static 
    double mean(const vector<Type>& values)
    {
        if (values.size() == 0)
            return 0;
        double mean = 0;
        for (size_t i = 0; i < values.size(); ++i)
            mean += values[i];
        mean /= values.size();
        return mean;
    }

    /*
     * Calculate the variance of a group of values
     */
    inline static 
    double var(const vector<Type>& values)
    {
        if (values.size() == 0)
            return 0;

        double mean = 0;
        for (size_t i = 0; i < values.size(); ++i)
            mean += values[i];
        mean /= values.size();

        double var = 0;
        for (size_t i = 0; i < values.size(); ++i)
            var += (values[i] - mean) * (values[i] - mean);
        var /= values.size() - 1;
        return var;
    }

    /*
     * Calculate the variance of a group of values
     */
    inline static 
    Type median(vector<Type>& values)
    {
        sort(values.begin(), values.end());
        return values[values.size() / 2];
    }

	inline static void 
	Switch (Type& A, Type& B)
	{
		Type temp = A;
		A = B;
		B = temp;
		return ;
	}		/* -----  end of method Utility::Switch  ----- */

	static void SortByIndex(vector<Type>& values, vector<int>& sortedIndex)
	{
		vector<bool> index_check;
		index_check.resize(sortedIndex.size());
		index_check.assign(sortedIndex.size(), false);
		for (size_t i = 0; i < sortedIndex.size(); i++)
			index_check[sortedIndex[i]] = true;
		for (size_t i = 0; i < sortedIndex.size(); i++)
			if (!index_check[sortedIndex[i]])
				cerr << __func__ << " ERROR invalid sortedIndex " << endl;

		if (sortedIndex.size() != values.size())
		{
			cerr << __func__ << " ERROR " << " : Two arrays have different length" << endl;
			return;
		}
		vector<Type> dups = values;
		for (size_t i = 0; i < sortedIndex.size(); i++)
			values[i] = dups[sortedIndex[i]];
	}
    
    // Equally partition a vector into specified number of groups
    // Take sum of individuals to form groups
    static void EqualGroupSum(const vector<Type>& values, vector<Type>& grouped_values, int group_cnt)
    {
        grouped_values.assign(group_cnt, 0);
        double each = (double)values.size() / group_cnt;
        for (size_t i = 0; i < values.size(); ++i)
        {
            if ((double)i >= each * group_cnt) break;
            grouped_values[(int)(i/each)] += values[i];
        }
    }

    static
    void Shuffle(vector<Type>& values)
    {
        const gsl_rng_type * T;
        gsl_rng * r;
        gsl_rng_env_setup();
        T = gsl_rng_default;
        r = gsl_rng_alloc (T);
        for (size_t i = 0; i < values.size(); ++i)
        {
            double u = gsl_rng_uniform (r);
            int idx = i + (int)(u * (values.size() - i));

            Type temp = values[i];
            values[i] = values[idx];
            values[idx] = temp;
        }
        gsl_rng_free (r);
    }

    /*
     * Values has been sorted from least to greatest.
     * If there is no item with values[i] = value,
     * return the i or i+1 such that values[i] < value < values[i+1]
     */
    static
    int BinarySearch(vector<Type>& values, const Type& value)
    {
        int low = 0;
        int high = values.size();
        int mid = 0;
        while (low < high)
        { 
            mid = (low + high) / 2;
            if (values[mid] > value)
                high = mid - 1;
            else if (values[mid] < value)
                low = mid + 1;
            else
                return mid;
        }
        mid = (low + high) / 2;
        return mid;
    }
};

template <typename Type, typename TComp = greater<Type> >
class UtilityTempComp
{
public:
	static void Sort(vector<Type>& keys, vector<int>& sortedIndex)
	{
		size_t i;

		sortedIndex.resize(keys.size());
		for (i = 0; i < sortedIndex.size(); i++)
			sortedIndex[i] = i;

		typedef multimap<Type, int, TComp> SORTMAP;
		typedef typename multimap<Type, int, TComp>::value_type PAIR;

		SORTMAP sortMap;

		for (i = 0; i < keys.size(); i++)
			sortMap.insert(PAIR(keys[i], i));

		i = 0;
		for_each_ele_in_group_temp(Iter, SORTMAP, sortMap)
		{
			keys[i] = Iter->first;
			sortedIndex[i++] = Iter->second;
		}
	}
};

#endif
