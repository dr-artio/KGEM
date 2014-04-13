package edu.gsu.cs.kgem.io

import edu.gsu.cs.kgem.model.Genotype
import java.io.{File, PrintStream}
import org.biojava3.core.sequence.DNASequence
import org.biojava3.core.sequence.io.FastaWriterHelper
import collection.JavaConversions._
import edu.gsu.cs.kgem.exec._

//import com.apporiented.algorithm.clustering._

import scala.Some

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:23 PM
 */
object OutputHandler {
  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   * { @see PrintStream} object, either file or stdout
   * @param gens
   * Collection of haplotypes (Result)
   */
  @deprecated
  def outputResult(out: PrintStream, gens: Iterable[Genotype]) = {
    outputHaplotypes(out, gens)
  }

  /**
   * Output corrected reads into specified {@see PrintStream}
   * @param out
   * { @see PrintStream} object, either file or stdout
   * @param gens
   * Collection of haplotypes (Result)
   * @param n
   * Number of reads
   */
  def outputResult(out: PrintStream, gens: Iterable[Genotype], n: Int, clean: (String => String) = (s => s)) = {
    val gg = gens.toIndexedSeq.sortBy(g => -g.freq)
    var i = 0
    val haplSeqs = gg.map(g => {
      val fn = (g.freq * n).asInstanceOf[Int]
      val cleanedSeq = trim(clean(g.toIntegralString), 'N')
      (0 until fn).map(x => {
        val dna = new DNASequence(cleanedSeq)
        dna.setOriginalHeader("read%d".format(i))
        i += 1
        dna
      })
    }).flatten
    writeFasta(out, haplSeqs)
  }

  /**
   * Output haplotypes into specified {@see PrintStream}
   * @param out
   * { @see PrintStream} object, either file or stdout
   * @param gens
   * Collection of haplotypes (Result)
   */
  def outputHaplotypes(out: PrintStream, gens: Iterable[Genotype], clean: (String => String) = (s => s)) = {
    val gg = gens.view.toIndexedSeq.sortBy(g => -g.freq)
    val haplSeqs = gg.map(g => {
      val seq = new DNASequence(trim(clean(g.toIntegralString), 'N'))
      seq.setOriginalHeader("haplotype%d_freq_%.10f".format(g.ID, g.freq))
      seq
    })
    writeFasta(out, haplSeqs)
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

  /**
   * Setup the output directory and files.
   * @param dir
   * The output directory file. This is where the reconstructed
   * haplotypes and corrected reads will be stored.
   * @return
   * Returns None if the output directory, or output files cannot
   * be created. If they are created successfully then it returns
   * Some((hapOutput: PrintStream, resultsOutput: PrintStream)).
   *
   */
  def setupOutput(dir: File): Option[(PrintStream)] = {
    // Try to make the output directory. If it fails, return None.
    val tmp = if (dir.getParentFile == null) new File(System.getProperty(USER_DIR)) else dir
    if (!tmp.exists()) {
      if (!tmp.mkdir()) {
        println("Cannot create output directory!")
        return None
      }
    }

    try {
      val out = new PrintStream(dir)
      return new Some[PrintStream](out)
    } catch {
      case _: Throwable => println("Cannot create file: " + dir.getAbsolutePath); return None
    }
  }

  /**
   * Merge two maps with collections as
   * values. Result map will contain union
   * of keys as new keyset and joint lists of
   * values
   * @param m1
   * Operand 1
   * @param m2
   * Operand 2
   * @tparam K
   * Generic parameter type of Key
   * @tparam V
   * Generic parameter type of Value
   * @return
   */
  def merge[K, V](m1: Map[K, Iterable[V]], m2: Map[K, Iterable[V]]): Map[K, Iterable[V]] = {
    val k1 = Set(m1.keysIterator.toList: _*)
    val k2 = Set(m2.keysIterator.toList: _*)
    val intersection = k1 & k2

    val r1 = for (key <- intersection) yield (key -> (m1(key) ++ m2(key)))
    val r2 = m1.filterKeys(!intersection.contains(_)) ++ m2.filterKeys(!intersection.contains(_))
    r2 ++ r1
  }

  def mean[T <% Double](vals: Iterable[T]) = {
    val N = vals.size
    vals.map(_.toDouble).sum / N
  }

  def sigma[T <% Double](vals: Iterable[T]) = {
    val meanValue = mean(vals)
    Math.sqrt(mean(vals.map(x => Math.pow(x - meanValue, 2))))
  }
}
