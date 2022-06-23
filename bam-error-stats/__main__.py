"""Report alignment error and concordance stats from a BAM file.

Helps with analyzing the output of CCS and DeepConsensus.
"""

from collections.abc import Sequence

from dataclasses import dataclass

from multiprocessing import Pool
import argparse
import math

import pysam


@dataclass
class AlnStats:
  read_name: str
  q_len: int
  effective_cov: float
  subread_passes: int
  pred_concordance: float
  supplementary: bool
  mapq: int
  concordance: float
  concordance_qv: float

  @staticmethod
  def header() -> str:
    return "#read,readLengthBp,effectiveCoverage,subreadPasses,predictedConcordance,alignmentType,alignmentMapq,hcReadLengthBp,concordance,concordanceQv,mismatchBp,nonHpInsertionBp,nonHpDeletionBp,hpInsertionBp,hpDeletionBp"

  def to_csv(self) -> str:
    aln_type = "Supplementary" if self.supplementary else "Primary"
    return f"{self.read_name},{self.q_len},{self.effective_cov:0.2f},{self.subread_passes},{self.pred_concordance:0.6f},{aln_type},{self.mapq},,{self.concordance:0.6f},{self.concordance_qv:0.2f},,,,,,,"


def concordance(record: pysam.AlignedSegment) -> float:
  stats, _ = record.get_cigar_stats()
  return float(stats[7]) / (stats[0] + stats[1] + stats[2])


def concordance_qv(concordance: float, read_len: int, has_errors: bool) -> float:
  qv_cap = 10.0 * math.log10(read_len + 1.0)
  if has_errors:
    return min(-10.0 * math.log10(1.0 - concordance), qv_cap)
  else:
    return qv_cap


def errors(record: pysam.AlignedSegment) -> int:
  stats, _ = record.get_cigar_stats()
  return stats[8] + stats[1] + stats[2]


def get_aln_stats(record: pysam.AlignedSegment) -> AlnStats or None:
  if record.is_unmapped() or not record.is_secondary():
    return None

  effective_cov = record.get_tag("ec") if record.has_tag("ec") else 0
  subread_passes = record.get_tag("np") if record.has_tag("np") else 0
  pred_concordance = record.get_tag("rq") if record.has_tag("rq") else 0
  concordance = concordance(record)
  concordance_qv = concordance_qv(concordance, record.infer_read_length(), errors(record) > 0)
  return AlnStats(record.query_name, record.query_length(), effective_cov, subread_passes, pred_concordance, record.is_supplementary(), record.mapping_quality(), concordance, concordance_qv)


def run(input_path: str, aln_stats_path: str, processes: int, chunk_size: int):
  with pysam.AlignmentFile(input_path, "rb", check_sq=False) as bam_file:
    with open(aln_stats_path, "w") as aln_stats_file:
      aln_stats_file.write(AlnStats.header())
      aln_stats_file.write("\n")

      with Pool(processes=processes) as p:
        for aln_stats in p.imap_unordered(get_aln_stats, bam_file.fetch(),
                                          chunk_size):
          if aln_stats is None:
            continue

          aln_stats_file.write(aln_stats.to_csv())
          aln_stats_file.write("\n")


if __name__ == "__main__":
  parser = argparse.ArgumentParser(
      description="Report alignment error and concordance stats from a BAM file.",
      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument("input", help="Input BAM file")
  parser.add_argument(
      "aln_stats", help="Output CSV file with per alignment statistics")
  parser.add_argument("-t", "--threads", type=int, default=4, help="Threads")
  parser.add_argument(
      "--chunk-size", type=int, default=1024, help="Chunk size per thread")
  args = parser.parse_args()
  run(args.input, args.aln_stats, args.threads, args.chunk_size)
