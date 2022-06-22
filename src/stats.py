"""Report alignment error and concordance stats from a BAM file.

Helps with analyzing the output of CCS and DeepConsensus.
"""

from collections.abc import Sequence

from multiprocessing import Pool
import argparse

import pysam


@dataclass
class AlnStats:
  length: int
  concordance: float

  def to_csv(self) -> str:
    return f"{self.length}, {self.concordance}"


def concordance(record: pysam.AlignedSegment) -> float:
  stats, _ = record.get_cigar_stats()
  return float(stats[7]) / (stats[0] + stats[1] + stats[2])


def get_aln_stats(record: pysam.AlignedSegment) -> AlnStats:
  return AlnStats(record.infer_read_length(), concordance(record))


def run(input_path: str, aln_stats_path: str, processes: int, chunk_size: int):
  with pysam.AlignmentFile(input_path, "rb") as bam_file:
    with open(aln_stats_path, "w") as aln_stats_file:
      with Pool(processes=processes) as p:
        for aln_stats in p.imap_unordered(get_aln_stats, bam_file.fetch(), chunk_size):
          aln_stats_file.write(aln_stats.to_csv())
          aln_stats_file.write("\n")


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Report alignment error and concordance stats from a BAM file.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("input", help="Input BAM file")
  parser.add_argument("aln_stats", help="Output CSV file with per alignment statistics")
  parser.add_argument("-t", "--threads", type=int, default=4, help="Threads")
  parser.add_argument("--chunk-size", type=int, default=1024, help="Chunk size per thread")
  args = parser.parse_args()
  run(args.input, args.aln_stats, args.threads, args.chunk_size)
