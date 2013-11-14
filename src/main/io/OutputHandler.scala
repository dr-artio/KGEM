package edu.gsu.cs.kgem.io

import edu.gsu.cs.kgem.model.{Genotype, Read}
import edu.gsu.cs.kgem.model.initialization.MaxDistanceSeedFinder.hammingDistance
import java.io.{File, PrintStream}
import org.biojava3.core.sequence.DNASequence
import org.biojava3.core.sequence.io.{FastaReaderHelper, FastaWriterHelper}
import collection.JavaConversions._
import org.apache.commons.lang3.StringUtils.getLevenshteinDistance

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:23 PM
 */
object OutputHandler {

  val UNCLUSTERED = -1

  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   */
  @deprecated
  def outputResult(out: PrintStream, gens: Iterable[Genotype]) = {
    outputHaplotypes(out, gens)
  }

  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   * @param n
   *          Number of reads
   */
  def outputResult(out: PrintStream, gens: Iterable[Genotype], n: Int, clean:(String => String) = (s => s)) = {
    val gg = gens.toIndexedSeq.sortBy(g => -g.freq)
    val haplSeqs = gg.map(g => {
      val fn = (g.freq * n).asInstanceOf[Int]
      val cleanedSeq = clean(g.toIntegralString)
      List.fill(fn)((cleanedSeq, g.freq))
    }).flatten.zipWithIndex.map(g => {
      val dna = new DNASequence(g._1._1)
      dna.setOriginalHeader("read%d_freq_%.10f".format(g._2,g._1._2))
      dna
    })
    writeFasta(out, haplSeqs)
  }

  /**
   * Output haplotypes into specified {@see PrintStream}
   * @param out
   *            {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of haplotypes (Result)
   */
  def outputHaplotypes(out: PrintStream, gens: Iterable[Genotype], clean:(String => String) = (s => s)) = {
    val gg = gens.toIndexedSeq.sortBy(g => -g.freq)
    val haplSeqs = gg.map(g => {
      val seq = new DNASequence(clean(g.toIntegralString))
      seq.setOriginalHeader("haplotype%d_freq_%.10f".format(g.ID,g.freq))
      seq
    })
    writeFasta(out, haplSeqs)
  }

  /**
   * Output reads with clustering info
   * @param out
   *             {@see PrintStream} object, either file or stdout
   * @param gens
   *             Collection of Haplotypes
   * @param reads
   *              Collection of reads
   * @param pqrs
   *             Clustering coefficients
   */
  def outputClusteredFasta(out: PrintStream, gens: Iterable[Genotype], reads: Iterable[Read],
                           pqrs: Array[Array[Double]], readsFile: File) = {
    val rs = FastaReaderHelper.readFastaDNASequence(readsFile)
    val genSeqs = gens.map(g => (g, g.toIntegralString.replaceAll("-",""))).toList

    println(gens.map(g => g.ID))

    val groups = reads.zipWithIndex.groupBy(rd => {
      val cl = column(pqrs, rd._2)
      val index = cl.indexWhere(p => p == cl.max)
      genSeqs(index)._1
    }).toMap

    val fReads = groups.map(rd => {
      rd._2.map(r => {
        val ds = rs(r._1.ids.head)
        ds.setOriginalHeader("h%d_%s_%d".format(rd._1.ID, ds.getOriginalHeader, r._1.ids.size))
        ds
      })
    }).flatten

    writeFasta(out, fReads)
  }

  private def trim(str: String, char: Char): String = {
    str.dropWhile(c => c == char).reverse.dropWhile(c => c == char).reverse
  }

  private def column[A, M[_]](matrix: M[M[A]], colIdx: Int)
                     (implicit v1: M[M[A]] => Seq[M[A]], v2: M[A] => Seq[A]): Seq[A] =
    matrix.map(_(colIdx))

  private def writeFasta(out: PrintStream, seq: Iterable[DNASequence]) {
    try {
      FastaWriterHelper.writeNucleotideSequence(out, seq)
    } finally {
      out.close()
    }
  }

  private def clusteringString(ggs: Iterable[(Genotype, Int)], pqrs: Array[Array[Double]], readIndex: Int) = {
    val sb = new StringBuffer()
    for (g <- ggs) {
      sb.append("_h%d=%s".format(g._1.ID, pqrs(g._2)(readIndex).toString))
    }
    sb.toString
  }

  /**
   * Setup the output directory and files.
   * @param dir
   *         The output directory file. This is where the reconstructed
   *         haplotypes and corrected reads will be stored.
   * @return
   *         Returns None if the output directory, or output files cannot
   *         be created. If they are created successfully then it returns
   *         Some((hapOutput: PrintStream, resultsOutput: PrintStream)).
   *
   */
  def setupOutputDir(dir: File): Option[(PrintStream, PrintStream, PrintStream, PrintStream, PrintStream)] = {
    // Try to make the output directory. If it fails, return None.
    if (!dir.exists()) {
      if (!dir.mkdir()) {
        println("Cannot create output directory!")
        return None
      }
    }

    // Try to open output files. If they fail, return None.
    val baseName = dir.getAbsolutePath() + File.separator
    val hapOutputName = "%s%s".format(baseName, "haplotypes.fas")
    val cleanedHapOutputName = "%s%s".format(baseName, "haplotypes_cleaned.fas")
    val readsOutputName = "%s%s".format(baseName, "reads.fas")
    val cleanedReadsOutputName = "%s%s".format(baseName, "reads_cleaned.fas")
    val readsClustered = "%s%s".format(baseName, "reads_clustered.fas")

    try {
        val hapout = new PrintStream(hapOutputName)
        try {
          val readsout = new PrintStream(readsOutputName)
          try {
            val hapclout = new PrintStream(cleanedHapOutputName)
            try {
              val readclsout = new PrintStream(cleanedReadsOutputName)
              try {
                val readsclust = new PrintStream(readsClustered)
                return Some((hapout, hapclout, readsout, readclsout, readsclust))
              } catch {
                case _: Throwable => println("Cannot create file: " + readsClustered); return None
              }
            } catch {
              case _: Throwable => println("Cannot create file: " + cleanedReadsOutputName); return None
            }
          } catch {
            case _: Throwable => println("Cannot create file: " + cleanedHapOutputName); return None
          }
        } catch {
          case _: Throwable => println("Cannot create file: " + readsOutputName); return None
        }
    } catch {
      case _: Throwable => println("Cannot create file: " + hapOutputName); return None
    }
  }

  /**
   * Merge two maps with collections as
   * values. Result map will contain union
   * of keys as new keyset and joint lists of
   * values
   * @param m1
   *           Operand 1
   * @param m2
   *           Operand 2
   * @tparam K
   *           Generic parameter type of Key
   * @tparam V
   *           Generic parameter type of Value
   * @return
   */
  def merge[K,V](m1: Map[K,Iterable[V]], m2: Map[K,Iterable[V]]): Map[K,Iterable[V]] = {
    val k1 = Set(m1.keysIterator.toList: _*)
    val k2 = Set(m2.keysIterator.toList: _*)
    val intersection = k1 & k2

    val r1 = for(key <- intersection) yield (key -> (m1(key) ++ m2(key)))
    val r2 = m1.filterKeys(!intersection.contains(_)) ++ m2.filterKeys(!intersection.contains(_))
    r2 ++ r1
  }
}
