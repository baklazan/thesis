#include<bits/stdc++.h>
using namespace std;

const double pi = 3.14159265358979323846;

struct Probability
{
  double log_val;
  
  Probability()
  {
    log_val = -numeric_limits<double>::infinity();
  }
  
  
  Probability(double logp)
  {
    log_val = logp;
  }
  
  static Probability fromP(double p)
  {
    return Probability(log(p));
  }
  
  double get()
  {
    return exp(log_val);
  }
  
  Probability& operator*=(const Probability& other)
  {
    log_val += other.log_val;
    return *this;
  }
  
  bool is_nan()
  {
    return log_val != log_val;
  }
};

Probability operator*(Probability a, Probability b)
{
  if(Probability(a.log_val + b.log_val).is_nan()) printf("%lf * %lf = nan\n", a.log_val, b.log_val);
  return Probability(a.log_val + b.log_val);
}

Probability operator+(Probability a, Probability b)
{
  if(a.log_val < b.log_val) swap(a, b);
  if(b.log_val == -numeric_limits<double>::infinity()) return a;
  return Probability(a.log_val + log(1 + exp(b.log_val - a.log_val)));
}

Probability operator/(Probability a, Probability b)
{
  return Probability(a.log_val - b.log_val);
}

bool operator<(Probability a, Probability b)
{
  return a.log_val < b.log_val;
}

Probability join(Probability a, Probability b)
{
  if((a + b).is_nan()) printf("%lf + %lf = nan\n", a.log_val, b.log_val);
  return a + b;
}

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
  
  int kmer_id_from_center(int *ref)
  {
    return kmer_id(&ref[-central_pos]);
  }
  
  Probability chance_of_signal(double v, int kmer)
  {
    double diff = v - mean[kmer];
    return Probability(add_const[kmer] - diff * diff * mult_const[kmer]);
  }
  
  Probability chance_mixture(double v, int kmer1, int kmer2)
  {
    return (chance_of_signal(v, kmer1) + chance_of_signal(v, kmer2)) / Probability::fromP(2);
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


Probability reresquiggle_flashback(double *signal, int signal_length, vector<int> kmer_ids, int start_buffer, int end_buffer, KmerModel *kmer_model, int min_event_length, double penalty)
{
  vector<vector<Probability> > dp(kmer_ids.size()*2, vector<Probability>(signal_length+1 , Probability::fromP(0)));
  for(int x=0; x<=start_buffer; x++)
  {
    dp[0][x] = Probability(-x * penalty);
  }
  
  for(int y=0; y<kmer_ids.size(); y++)
  {
    for(int x=min_event_length; x<=signal_length; x++)
    {
      Probability start_price(0);
      for(int i=1; i<=min_event_length; i++)
      {
        start_price *= kmer_model->chance_of_signal(signal[x-i], kmer_ids[y]);
      }
      Probability continue_price = kmer_model->chance_of_signal(signal[x-1], kmer_ids[y]);
      dp[y*2+1][x] = join(dp[y*2+1][x-1] * continue_price, dp[y*2][x-min_event_length] * start_price);
    }
    if(y+1 < kmer_ids.size())
    {
      for(int x=1; x<=signal_length; x++)
      {
        Probability continue_price = kmer_model->chance_mixture(signal[x-1], kmer_ids[y], kmer_ids[y+1]);
        dp[y*2+2][x] = join(dp[y*2+1][x], dp[y*2+2][x-1] * continue_price);
      }
    }
  }
  Probability res = Probability::fromP(0);
  for(int x=signal_length - end_buffer; x <= signal_length; x++)
  {
    res = join(res, dp.back()[x] / Probability((signal_length-x)*penalty));
  }
  return res;
}

Probability reresquiggle(double *signal, int signal_length, vector<int> kmer_ids, int start_buffer, int end_buffer, KmerModel *kmer_model, int min_event_length, double penalty)
{
  vector<vector<Probability> > dp(kmer_ids.size()+1, vector<Probability>(signal_length+1 , Probability::fromP(0)));
  for(int x=0; x<=start_buffer; x++)
  {
    dp[0][x] = Probability(-x * penalty);
  }
  
  for(int y=0; y<kmer_ids.size(); y++)
  {
    for(int x=min_event_length; x<=signal_length; x++)
    {
      Probability start_price(0);
      for(int i=1; i<=min_event_length; i++)
      {
        start_price *= kmer_model->chance_of_signal(signal[x-i], kmer_ids[y]);
      }
      Probability continue_price = kmer_model->chance_of_signal(signal[x-1], kmer_ids[y]);
      dp[y+1][x] = join(dp[y+1][x-1] * continue_price, dp[y][x-min_event_length] * start_price);
    }
  }
  Probability res = Probability::fromP(0);
  for(int x=signal_length - end_buffer; x <= signal_length; x++)
  {
    res = join(res, dp.back()[x] / Probability((signal_length-x)*penalty));
  }
  return res;
}

class full_computation
{
  KmerModel *kmer_model;
  ResquiggledRead *read;
  int half_bandwidth, *reference, reference_length, patch_size;
  double *result;
  int height;
  vector<int> row_start, row_end, kmer;
  bool flashbacks;
  int M;
  Probability SNP_prior, nonSNP_prior;
  
  
public:
  full_computation(int *_reference, int reference_length, KmerModel *kmer_model, ResquiggledRead *read, 
                   int half_bandwidth, double *result, int M, bool flashbacks, double expected_SNPs)
  {
    this->kmer_model = kmer_model;
    this->read = read;
    this->half_bandwidth = half_bandwidth;
    this->result = result;
    reference = _reference + read->start_in_reference;
    this->reference_length = reference_length;
    this->flashbacks = flashbacks;
    height = read->event_start.size();
    this->M = M;
    
    if(flashbacks)
    {
      patch_size = kmer_model->k + 1;
    }
    else
    {
      patch_size = kmer_model->k;
    }
    
    SNP_prior =  Probability::fromP(expected_SNPs/3);
    nonSNP_prior =  Probability::fromP(1-expected_SNPs);
    
    row_start = vector<int>(height+1);
    row_end = vector<int>(height+1);
    row_start[0] = read->event_start[0];
    row_end[0] = read->event_start[0] + 1;
    for(int i=0; i<height; i++)
    {
      row_start[i+1] = read->event_start[max(0, i-half_bandwidth)] + 1;
      row_end[i+1] = read->event_end[min(height-1, i+half_bandwidth)] + 1;
    }
    kmer = vector<int> (height);
    for(int i=0; i<height; i++)
    {
      kmer[i] = kmer_model->kmer_id_from_center(&(reference[i]));
    }
  }
  
  void compute_flashback_layer_backward(int y, Probability *layer, Probability *row, int previous_kmer, int current_kmer)
  {
    int st = row_start[y];
    int x = row_end[y] - 1;
    layer[x - st] = row[x - st];
    for(int x = row_end[y] - 2; x >= st; x--)
    {
      layer[x - st] = join(row[x - st], layer[x + 1 - st] * kmer_model->chance_mixture(read->signal[x], previous_kmer, current_kmer));
    }
  }
  
  Probability cost_M_steps(int signal_start_id, int kmer_id)
  {
    Probability res(0);
    for(int i=0; i<M; i++)
    {
      res *= kmer_model->chance_of_signal(read->signal[signal_start_id + i], kmer_id);
    }
    if(res.is_nan()) printf("product is nan\n");
    return res;
  }
  
  void compute_row_backward(int y, Probability *current_row, Probability *next_row, int current_kmer, int next_kmer = -1)
  {
    vector<Probability> inter_layer;
    if(next_row != nullptr && flashbacks)
    {
      inter_layer.resize(row_end[y+1] - row_start[y+1]);
      compute_flashback_layer_backward(y+1, &inter_layer[0], next_row, current_kmer, next_kmer);
      next_row = &inter_layer[0];
    }
    
    int st = row_start[y];
    int x = row_end[y] - 1;
    if(next_row == nullptr)
    {
      current_row[x - st] = Probability(0);
    }
    else
    {
      if(x + M >= row_start[y+1] && x+M < row_end[y+1])
      {
        current_row[x - st] = next_row[x+M - row_start[y+1]] * cost_M_steps(x, next_kmer);
      }
      else
      {
        current_row[x - st] = Probability::fromP(0);
      }
    }
    for(int x = row_end[y] - 2; x >= st; x--)
    {
      current_row[x - st] = current_row[x+1 - st] * kmer_model->chance_of_signal(read->signal[x], current_kmer);
      if(next_row != nullptr && x+M >= row_start[y+1] && x+M < row_end[y+1])
      {
        Probability alternative = next_row[x+M-row_start[y+1]] * cost_M_steps(x, next_kmer);
        current_row[x-st] = join(current_row[x-st], alternative);
      }
    }
  }
  
  void compute_flashback_layer_forward(int y, Probability *layer, Probability *row, int current_kmer, int next_kmer)
  {
    int st = row_start[y];
    row[0] = layer[0];
    for(int x=st+1; x<row_end[y]; x++)
    {
      row[x - st] = join(layer[x - st], row[x - 1 - st] * kmer_model->chance_mixture(read->signal[x-1], current_kmer, next_kmer));
    }
  }
  
  void compute_row_forward(int y, Probability *current_row, Probability *previous_row, int current_kmer, int next_kmer = -1)
  {
    Probability *row = current_row;
    vector<Probability> inter_layer;
    if(flashbacks && next_kmer != -1)
    {
      inter_layer.resize(row_end[y] - row_start[y]);
      row = &inter_layer[0];
    }
    int st = row_start[y];
    for(int x=st; x<row_end[y]; x++)
    {
      row[x - st] = Probability::fromP(0);
      if(x > st)
      {
        if(kmer_model->chance_of_signal(read->signal[x-1], current_kmer).is_nan()) printf("nans, nans everywhere\n");
        row[x - st] = row[x - 1 - st] * kmer_model->chance_of_signal(read->signal[x-1], current_kmer);
        if(row[x - st].is_nan()) printf("producing nan A\n");
      }
      if(x-M >= row_start[y-1] && x-M < row_end[y-1])
      {
        if(row[x - st].is_nan()) printf("already nan B\n");
        row[x - st] = join(row[x - st], previous_row[x-M - row_start[y-1]] * cost_M_steps(x-M, current_kmer));
        if(row[x - st].is_nan()) printf("producing nan B\n");
      }
    }
    for(int x=st; x<row_end[y]; x++) if(row[x - st].is_nan()) printf("nan at intermediate row %d, column %d before flashback layer\n", x, y);
    if(flashbacks && next_kmer != -1)
    {
      compute_flashback_layer_forward(y, row, current_row, current_kmer, next_kmer);
    }
    for(int x=st; x<row_end[y]; x++) if(current_row[x - st].is_nan()) printf("nan at row %d, column %d\n", x, y);
    for(int x=st; x<row_end[y]; x++) if(row[x - st].is_nan()) printf("nan at intermediate row %d, column %d\n", x, y);
  }
  
  void compute()
  {
    for(int i=0; i<height*4; i++)
    {
      result[i] = -numeric_limits<double>::infinity();
    }
    
    vector<vector<Probability> > cost_to_end(height+1);
    cost_to_end.back().resize(row_end[height] - row_start[height]);
    compute_row_backward(height, &cost_to_end.back()[0], nullptr, kmer[height-1]);
    for(int y=height-1; y>0; y--)
    {
      vector<vector<Probability> > candidates(1, vector<Probability>(row_end[y] - row_start[y]));
      compute_row_backward(y, &candidates[0][0], &cost_to_end[y+1][0], kmer[y-1], kmer[y]);
      if(y + patch_size <= height)
      {
        int position = y + patch_size - 1 - kmer_model->central_pos;
        int original = reference[position];
        for(int i=0; i<4; i++)
        {
          if(i == original) continue;
          reference[position] = i;
          vector<vector<Probability> > patch(patch_size);
          for(int yy = y+patch_size-1; yy >= y; yy--)
          {
            patch[yy-y].resize(row_end[yy] - row_start[yy]);
            Probability *next = yy+1-y < patch_size? &patch[yy+1-y][0] : &cost_to_end[yy+1][0];
            int current_kmer = kmer_model->kmer_id_from_center(&reference[yy-1]);
            int next_kmer = kmer_model->kmer_id_from_center(&reference[yy]);
            compute_row_backward(yy, &patch[yy-y][0], next, current_kmer, next_kmer);
          }
          candidates.push_back(patch[0]);
        }
        reference[position] = original;
      }
      
      
      cost_to_end[y] = vector<Probability>(row_end[y] - row_start[y], Probability::fromP(0));
      for(int i=0; i<1; i++) //candidates.size(); i++)
      {
        for(int j=0; j<cost_to_end[y].size(); j++)
        {
          cost_to_end[y][j] = max(cost_to_end[y][j], candidates[i][j]);
        }
      }
    }
    

    vector<vector<Probability> > cost_from_start(height+1);
    cost_from_start[0] = vector<Probability>(1, 0);
    for(int y=1; y<=height; y++)
    {
      vector<vector<Probability> > candidates(1, vector<Probability>(row_end[y] - row_start[y]));
      compute_row_forward(y, &candidates[0][0], &cost_from_start[y-1][0], kmer[y-1], y < height? kmer[y] : -1);
      if(y >= patch_size)
      {
        int position = y - 1 - kmer_model->central_pos;
        int original = reference[position];
        for(int i=0; i<4; i++)
        {
          if(i == original) continue;
          reference[position] = i;
          Probability *previous = &cost_from_start[y-patch_size][0];
          vector<vector<Probability> > patch;
          for(int yy = y-patch_size+1; yy <= y; yy++)
          {
            patch.push_back(vector<Probability>(row_end[yy] - row_start[yy]));
            int current_kmer = kmer_model->kmer_id_from_center(&reference[yy-1]);
            int next_kmer = -1;
            if(yy < height) next_kmer = kmer_model->kmer_id_from_center(&reference[yy]);
            compute_row_forward(yy, &patch.back()[0], previous, current_kmer, next_kmer);
            previous = &patch.back()[0];
          }
          candidates.push_back(patch.back());
          Probability best = Probability::fromP(0);
          for(int x=row_start[y]; x<row_end[y]; x++)
          {
            best = join(best, patch.back()[x - row_start[y]] * cost_to_end[y][x - row_start[y]]);
          }
          best.log_val /= 10;
          result[position * 4 + i] = (best * SNP_prior).log_val;
          for(int yy= max(0, position - kmer_model->k + 1); yy<min(position + kmer_model->k, reference_length); yy++)
          {
            if(yy == position) continue;
            if(best.is_nan()) printf("b\n");
            result[yy * 4 + reference[yy]] = (Probability(result[yy * 4 + reference[yy]]) + best * SNP_prior).log_val;
          }
        }
        reference[position] = original;
        Probability best = Probability::fromP(0);
        for(int x=row_start[y]; x<row_end[y]; x++)
        {
          if(candidates[0][x - row_start[y]].is_nan()) printf("d\n");
          if(cost_to_end[y][x - row_start[y]].is_nan()) printf("e\n");
          best = join(best, candidates[0][x - row_start[y]] * cost_to_end[y][x - row_start[y]]);
        }
        best.log_val /= 10;
        result[position * 4 + reference[position]] = (Probability(result[position * 4 + reference[position]]) + best * nonSNP_prior).log_val;
      }
      
      cost_from_start[y] = vector<Probability>(row_end[y] - row_start[y], Probability::fromP(0));
      for(int i=0; i<1; i++) //candidates.size(); i++)
      {
        for(int j=0; j<cost_from_start[y].size(); j++)
        {
          cost_from_start[y][j] = max(cost_from_start[y][j], candidates[i][j]);
        }
      }
    }
    

    
    for(int i=0; i<height; i++)
    {
      bool ok = false;
      for(int j=0; j<4; j++)
      {
        if(result[i*4+j] != -numeric_limits<double>::infinity())
        {
          ok = true;
        }
      }
      if(!ok)
      {
        for(int j=0; j<4; j++)
        {
          result[i*4+j] = 0;
        }
      }
    }
  }
  
};


extern "C"
{
  void compute_probabilities(int *reference, int reference_length, KmerModel *kmer_model, 
                             ResquiggledRead *read, int *interesting_positions, int interesting_count, double *result,
                             int min_event_length, int window_before, int window_after, int buffer_size, double penalty,
                             bool flashbacks, double expected_SNPs)
  {
    Probability (*reresquiggle_func)(double*, int, vector<int>, int, int, KmerModel*, int, double) = flashbacks ? reresquiggle_flashback : reresquiggle;
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
        double prior = letter == reference_value ? 1-expected_SNPs : expected_SNPs / 3;
        result[i * 4 + letter] = reresquiggle_func(&(read->signal[st_in_signal]),
                                                   en_in_signal - st_in_signal, 
                                                   kmer_ids, 
                                                   start_buffer, 
                                                   end_buffer, 
                                                   kmer_model, 
                                                   min_event_length,
                                                   penalty).log_val / 10 + log(prior);
      }
      reference[pos] = reference_value;
    }
  }
  
  void compute_all_probabilities(int *reference, int reference_length, KmerModel *kmer_model, 
                                 ResquiggledRead *read, int half_bandwidth, double *result, 
                                 int min_event_length, bool flashbacks, double expected_SNPs)
  {
    full_computation(reference, reference_length, kmer_model, read, half_bandwidth, result, min_event_length, flashbacks, expected_SNPs).compute();
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