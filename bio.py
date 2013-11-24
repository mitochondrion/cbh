from collections import defaultdict

def frequency_per_kmer(sequence, k):
  kmer_counts = defaultdict(int)

  for start_idx in range((len(sequence) - k) + 1):
    kmer_counts[sequence[start_idx: start_idx + k]] += 1

  return kmer_counts

###

def kmers_per_frequeny(sequence, k):
  count_kmers = defaultdict(list)
  kmer_counts = frequency_per_kmer(sequence, k)
  
  for kmer,count in kmer_counts.items():
    count_kmers[count].append(kmer)

  return count_kmers

###

def most_frequent_kmers(sequence, k):
  count_kmers = kmers_per_frequeny(sequence, k)
 
  return count_kmers[max(count_kmers.keys())]

###

def reverse_complement(sequence):
  complement_map = { 'A': 'T', 'T': 'A', 'G':'C', 'C': 'G' }
  complement = [complement_map[n] for n in sequence]
  complement.reverse()
  return "".join(complement)

###

def occurrences(pattern, sequence):
  occurrences = []
  pattern_len = len(pattern)

  for start_idx in range((len(sequence) - pattern_len) + 1):
    if pattern == sequence[start_idx: start_idx + pattern_len]:
      occurrences.append(start_idx)

  return occurrences

###

def clumping_kmers(sequence, k, window_size, min_occurrence):
  clumping_kmers = set()
  
  for start_idx in range((len(sequence) - window_size) + 1):
    window = sequence[start_idx: start_idx + window_size]
    kmer_counts = frequency_per_kmer(window, k)
    for kmer, count in kmer_counts.items():
      if count >= min_occurrence:
        clumping_kmers.add(kmer)

  return clumping_kmers

###

def faster_clumping_kmers(sequence, k, window_size, min_occurrence):
  sequence_len = len(sequence)
  clumping_kmers = set()
  window = sequence[0 : window_size]
  kmer_counts = frequency_per_kmer(window, k)

  for kmer, count in kmer_counts.items():
    if count >= min_occurrence:
      clumping_kmers.add(kmer)

  print("=====:", range(window_size + 1 - k, (sequence_len - k) + 1))

  for next_kmer_idx in range(window_size + 1 - k, (sequence_len - k) + 1):
    if next_kmer_idx % 100000 == 0:
      print(next_kmer_idx)
    
    previous_kmer_idx = next_kmer_idx - (window_size - k) - 1
    previous_kmer = sequence[previous_kmer_idx : previous_kmer_idx + k]
    next_kmer = sequence[next_kmer_idx : next_kmer_idx + k]
    
    kmer_counts[previous_kmer] -= 1
    kmer_counts[next_kmer] += 1

    for kmer in [previous_kmer, next_kmer]: 
      if kmer_counts[kmer] >= min_occurrence:
        clumping_kmers.add(kmer)

  return clumping_kmers

###

def gc_skews(genome):
  skew_increment = defaultdict(int, {'g': 1, 'G': 1, 'c': -1, 'C': -1})

  skews = [0]
  
  for base in genome:
    skews.append(skews[-1] + skew_increment[base])

  return skews

###

def min_gc_skews(genome):
  skew_increment = defaultdict(int, {'g': 1, 'G': 1, 'c': -1, 'C': -1})

  skew = 0
  min_skew = 0
  min_skews = []
  
  for idx, base in enumerate(genome, 1):
    skew = skew + skew_increment[base]
    if skew == min_skew:
      min_skews.append(idx)
    elif skew < min_skew:
      min_skew = skew
      min_skews = [idx]

  return min_skews
