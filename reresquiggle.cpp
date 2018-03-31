#include<bits/stdc++.h>
using namespace std;

const double pi = 3.14159265358979323846;

class KmerModel
{
	vector<double> add_const, mult_const, mean;
public:
	int k, central_pos;
	KmerModel(int _k, int _central_pos, double *_mean, double *sigma)
	{
		k = _k;
		central_pos = _central_pos;
		int kmer_count = 1;
		for(int i=0; i<k; i++) kmer_count *= 4;
		add_const.resize(kmer_count);
		mult_const.resize(kmer_count);
		mean.resize(kmer_count);
		for(int i=0; i<kmer_count; i++)
		{
			mean[i] = _mean[i];
			add_const[i] = log(1 / sqrt(2 * pi * sigma[i] * sigma[i]));
			mult_const[i] = 1 / (2 * sigma[i] * sigma[i]);
		}
	}
	
	int kmer_id(int *ref)
	{
		int res = 0;
		for(int i=0; i<k; i++)
		{
			res *= 4;
			res += ref[i];
		}
		return res;
	}
	
	double log_chance_of_signal(double v, int kmer)
	{
		double diff = v - mean[kmer];
		return add_const[kmer] - diff * diff * mult_const[kmer];
	}
};

struct ResquiggledRead
{
	vector<double> signal;
	vector<int> event_start, event_end;
	int start_in_reference;
	
	ResquiggledRead(double *_signal, int signal_length, int *_event_start, int *event_length, int event_count, int _start_in_reference)
	{
		start_in_reference = _start_in_reference;
		signal.resize(signal_length);
		for(int i=0; i<signal_length; i++)
		{
			signal[i] = _signal[i];
		}
		event_start.resize(event_count);
		event_end.resize(event_count);
		for(int i=0; i<event_count; i++)
		{
			event_start[i] = _event_start[i];
			event_end[i] = event_start[i] + event_length[i];
		}
	}
};

double log_average(double a, double b)
{
	if(a < b) swap(a, b);
	return a + log(1 + exp(b - a)) - log(2);
}

double reresquiggle(double *signal, int signal_length, vector<int> kmer_ids, int start_buffer, int end_buffer, KmerModel *kmer_model, int min_event_length, double penalty)
{
	vector<vector<double> > dp(kmer_ids.size()*2, vector<double>(signal_length+1 ,-numeric_limits<double>::infinity()));
	for(int x=0; x<start_buffer; x++)
	{
		dp[0][x] = -x * penalty;
	}
	
	for(int y=0; y<kmer_ids.size(); y++)
	{
		for(int x=min_event_length; x<=signal_length; x++)
		{
			double start_price = 0;
			for(int i=1; i<=min_event_length; i++)
			{
				start_price += kmer_model->log_chance_of_signal(signal[x-i], kmer_ids[y]);
			}
			double continue_price = kmer_model->log_chance_of_signal(signal[x-1], kmer_ids[y]);
			dp[y*2+1][x] = max(dp[y*2+1][x-1] + continue_price, dp[y*2][x-min_event_length] + start_price);
		}
		if(y+1 < kmer_ids.size())
		{
			for(int x=1; x<=signal_length; x++)
			{
				double continue_price = log_average(kmer_model->log_chance_of_signal(signal[x-1], kmer_ids[y]),
				                                    kmer_model->log_chance_of_signal(signal[x-1], kmer_ids[y+1]));
				dp[y*2+2][x] = max(dp[y*2+1][x], dp[y*2+2][x-1] + continue_price);
			}
		}
	}
	double res = -numeric_limits<double>::infinity();
	for(int x=signal_length - end_buffer; x <= signal_length; x++)
	{
		res = max(res, dp.back()[x] - (signal_length-x)*penalty);
	}
	return res;
}

extern "C"
{
	void compute_probabilities(int *reference, int reference_length, KmerModel *kmer_model, 
	                           ResquiggledRead *read, int *interesting_positions, int interesting_count, double *result,
	                           int min_event_length, int window_before, int window_after, int buffer_size, double penalty)
	{
		for(int i=0; i<interesting_count; i++)
		{
			int pos = interesting_positions[i];
			int st_in_ref = pos - window_before;
			int en_in_ref = pos + 1 + window_after;
			int context_st_in_ref = st_in_ref - kmer_model->central_pos;
			int context_en_in_ref = en_in_ref + kmer_model->k - kmer_model->central_pos - 1;
			if(context_st_in_ref < 0 || context_en_in_ref >= reference_length) continue;
			
			int st_in_read = st_in_ref - read->start_in_reference;
			int en_in_read = en_in_ref - read->start_in_reference;
			if(st_in_read - buffer_size < 0 || en_in_read + buffer_size > read->event_start.size()) continue;
			
			int st_in_signal = read->event_start[st_in_read - buffer_size];
			int start_buffer = read->event_start[st_in_read + buffer_size] - st_in_signal;
			int en_in_signal = read->event_end[en_in_read + buffer_size - 1];
			int end_buffer = en_in_signal - read->event_end[en_in_read - buffer_size - 1];
			
			int reference_value = reference[pos];
			for(int letter = 0; letter < 4; letter++)
			{
				reference[pos] = letter;
				vector<int> kmer_ids;
				for(int j=context_st_in_ref; j+kmer_model->k <= context_en_in_ref; j++)
				{
					kmer_ids.push_back(kmer_model->kmer_id(&reference[j]));
				}
				result[i * 4 + letter] = reresquiggle(&(read->signal[st_in_signal]),
				                                 en_in_signal - st_in_signal, 
				                                 kmer_ids, 
				                                 start_buffer, 
				                                 end_buffer, 
				                                 kmer_model, 
				                                 min_event_length,
				                                 penalty);
			}
			reference[pos] = reference_value;
		}
	}
	
	KmerModel *new_KmerModel(int k, int central_pos, double *mean, double *sigma)
	{
		return new KmerModel(k, central_pos, mean, sigma);
	}
	
	void delete_KmerModel(KmerModel *kmer_model)
	{
		delete kmer_model;
	}
	
	ResquiggledRead *new_ResquiggledRead(double *signal, int signal_length, int *event_start, int *event_length, int event_count, int start_in_reference)
	{
		return new ResquiggledRead(signal, signal_length, event_start, event_length, event_count, start_in_reference);
	}
	
	void delete_ResquiggledRead(ResquiggledRead *read)
	{
		delete read;
	}
}