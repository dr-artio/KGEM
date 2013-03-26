package edu.gsu.cs.kgem.exec

import collection.mutable
import edu.gsu.cs.kgem.model.{KGEM, Genotype, Read}
import net.sf.samtools.{SAMFileHeader, SAMRecord}
import scala.io.Source._
import scala.util.Random
import net.sourceforge.argparse4j.ArgumentParsers
import java.io._
import net.sourceforge.argparse4j.inf.ArgumentParserException

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/15/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
object Main {
  val READS_PARAMETER = "reads"
  val K_PARAMETER = "k"
  val OUTPUT_PARAMETER = "out"

  def main(args: Array[String]) {
    val (k, fl, out) = argumentParsing(args)

    val s = System.currentTimeMillis

    val lines = fromFile(fl).getLines
    val readsMap = toCounterMap(lines)
    val SAMRecords = toSAMRecords(readsMap.keys)
    val reads = toReads(SAMRecords)
    initReadFreqs(reads, readsMap)
    val seeds = sample(reads, k)
    var gens = seeds.map(s => new Genotype(s.seq)).toList
    KGEM.initReads(reads.toList)
    gens = KGEM.run(gens)
    for (g <- gens) {
      out.println(g.freq + "\n" + g.toIntegralString.replace("-", ""))
    }
    println("The whole procedure took " +
      ((System.currentTimeMillis - s) * 0.0001 / 6) + " minutes\n" +
      "Total number of haplotypes is " + gens.size + "\nbye bye")
  }

  /**
   * Method for choosing random sample of size @size from the
   * collection of objects. Returns the whole list is @size
   * is greater than size of the original collection.
   * @param iter
   *             Any iterable collection
   * @param size
   *             Size of the required sample
   * @tparam T
   *           Generic parameter (Class of objects in list)
   * @return
   *         List of randomly chosen elements from collection
   *         of specified size @size
   */
  def sample[T](iter: Iterable[T], size: Int) = {
    if (iter.size < size) iter.toList
    var res = new mutable.MutableList[T]()
    val rnd = new Random(System.currentTimeMillis)
    var needed = size
    var len = iter.size
    val iterator = iter.iterator
    while (needed > 0 && iterator.hasNext) {
      val item = iterator.next
      if (rnd.nextInt(len) < needed) {
        res += item
        needed -= 1
      }
      len -= 1
    }
    res
  }

  def initReadFreqs(reads: Iterable[Read], readsMap: mutable.Map[String, Int]) {
    reads foreach (r => {
      r.freq = readsMap(r.seq)
    })
  }

  def toReads(sams: Iterable[SAMRecord]) = {
    for (sam <- sams) yield new Read(sam)
  }

  def toCounterMap(lines: Iterator[String]): mutable.Map[String, Int] = {
    val readsMap = mutable.Map[String, Int]()
    for (l <- lines) {
      if (readsMap contains l) {
        readsMap(l) += 1
      } else {
        readsMap.put(l, 1)
      }
    }
    readsMap
  }

  def toSAMRecords(reads: Iterable[String]) = {
    for (r <- reads)
    yield toSAMRecord(r)
  }

  def toSAMRecord(st: String) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(0)
    sam.setReadString(st)
    sam
  }

  /**
   * Method for parsing command line parameters
   * and handle mistakes in parameters list
   * @param args
   *             Array of command line parameters
   * @return
   *         Tuple of parsed and converted parameters
   */
  private def argumentParsing(args: Array[String]) = {
    val parser = ArgumentParsers.newArgumentParser("KGEM")
      .description("Error correction based on KGEM.");

    var k = 50
    var fl: File = null
    var out = System.out

    parser.addArgument(READS_PARAMETER)
      .metavar("ReadsFile")
      .help("File containing preprocessed sequencing data"
      + " file with extension (.txt) "
      + "reads in extended format")
      .`type`(classOf[File])

    parser.addArgument("-k").dest(K_PARAMETER)
      .metavar("K")
      .setDefault[Integer](k)
      .`type`(classOf[Integer])
      .help("Parameter k - the size of sample being randomly chosen "
      + "as seeds. Depends on expectation of variability expected "
      + "on exploring region Default: " + k + ")");

    parser.addArgument("-o", "-out").dest(OUTPUT_PARAMETER)
      .metavar("O")
      .setDefault[PrintStream](out)
      .`type`(classOf[FileOutputStream])
      .help("Output file name. Default: system.out")

    try {
      val n = parser.parseArgs(args)
      k = n.getInt(K_PARAMETER)
      fl = n.get(READS_PARAMETER).asInstanceOf[File]
      val outO = n.get(OUTPUT_PARAMETER)
      if (outO.isInstanceOf[FileOutputStream]) out = new PrintStream(outO.asInstanceOf[OutputStream])
    } catch {
      case e: ArgumentParserException => {
        parser.handleError(e)
        System.exit(1)
      }
    }
    (k, fl, out)
  }
}
