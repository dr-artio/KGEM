package edu.gsu.cs.kgem

import edu.gsu.cs.kgem.model._
import collection.mutable
import edu.gsu.cs.kgem.io.{ArgumentParser, SAMParser}
import java.io.{PrintStream, File}
import net.sf.samtools.{SAMFileHeader, SAMRecord}
import org.biojava3.core.sequence.io.FastaReaderHelper._
import collection.JavaConversions._
import scala.io.Source.fromFile
import java.util.Date
import java.text.SimpleDateFormat
import edu.gsu.cs.kgem.io.OutputHandler._
import edu.gsu.cs.kgem.io.Config
import scala.Some
import org.biojava3.core.sequence.DNASequence
import edu.gsu.cs.kgem.io.Config
import scala.Some
import edu.gsu.cs.kgem.io.Config
import scala.Some

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:02 PM
 */
package object exec {
  private val sdf = new SimpleDateFormat("[hh:mm:ss a]")
  private val FASTA = Array[String](".fas", ".fa", ".fasta")
  val USER_DIR = "user.dir"
  val KGEM_STR = "kGEM version %s: Local Reconstruction for Mixed Viral Populations."
  val LINE = "------------------------------------------------------------"
  var hap: PrintStream = null
  var hapcl: PrintStream = null
  var res: PrintStream = null
  var rescl: PrintStream = null
  var clust: PrintStream = null
  var config: Config = null
  var reads: List[Read] = null
  var seqs: Map[String, DNASequence] = null
  var k: Int = -1
  var threshold: Int = 0
  var n: Int = 0

  protected[exec] def parseArgs(args: Array[String]) = {
    ArgumentParser.parseArguments(args) match {
      case None => sys.exit(1)
      case Some(config: Config) => {
        this.config = config
        setupOutputDir(config.output) match {
          case None => sys.exit(1)
          case Some((hap: PrintStream, hapcl: PrintStream, res: PrintStream, rescl: PrintStream, clust: PrintStream)) => {
            this.hap = hap
            this.hapcl = hapcl
            this.res = res
            this.rescl = rescl
            this.clust = clust
          }
        }
      }
    }
  }

  def initInputData(readsFile: File = config.readsFile) = {
    if (!FASTA.exists(readsFile.getName.toLowerCase.endsWith)) sys.exit(1)
    reads = initFastaReads(readsFile).toList
    k = if (config.clustering != null) config.k * 2 else config.k
    threshold = config.threshold
    n = reads.map(r => r.freq).sum.toInt
  }

  def executeKgem(reads: List[Read] = reads, k: Int = k, threshold: Int = threshold) = {
    KGEM.initReads(reads.toList)
    val gens = if (config.consensusFile == null) {
      if (config.prThr >= 0) KGEM.initThreshold(config.prThr)
      else KGEM.initThreshold
      val seeds = MaxDistanceSeedFinder.findSeeds(reads, k, threshold)
      KGEM.run(seeds)
    } else {
      val seeds = initFastaReads(config.consensusFile).map(r => new Genotype(r.seq))
      KGEM.run(seeds)
    }
    gens.toList
  }

  protected[exec] def clusteringHandler(gens: List[Genotype]) = {
    if (config.clustering != null) {
      val readSeqs = readFastaDNASequence(config.clustering).toMap
      val clusteredReads = KGEM.clustering(gens, readSeqs, config.k)
      writeFasta(clust, clusteredReads)
    }
  }

  protected[exec] def outputResults(gens: List[Genotype], s: Long) = {
    outputHaplotypes(hap, gens)
    outputHaplotypes(hapcl, gens, s => s.replaceAll("-", ""))
    outputResult(res, gens, n)
    outputResult(rescl, gens, n, s => s.replaceAll("-", ""))

    log("The whole procedure took %.2f minutes".format(((System.currentTimeMillis - s) * 0.0001 / 6)))
    log("Total number of haplotypes is %d".format(gens.size))
    log("bye bye")
  }


  /**
   * Read SAM file and fold according to extended
   * sequences.
   * @param fl
   * SAM file
   * @return
   * Iterable collection of Reads
   */
  @deprecated
  def initSAMReads(fl: File): Iterable[Read] = {
    val samRecords = SAMParser.readSAMFile(fl)
    var extSAMRecords = samRecords.map(s => SAMParser.toExtendedString(s))
    val l = extSAMRecords.map(s => s.length).max
    extSAMRecords = samRecords.map(s => SAMParser.toExtendedString(s, l))
    val samMap = samRecords.map(s => (SAMParser.toExtendedString(s, l), s)).toMap
    val readsMap = flip(samMap.map(s => (s._2.getReadName, s._1)).toMap)
    val reads = toReads(samMap)
    initReadFreqs(reads, readsMap)
    reads
  }

  /**
   * Old method for parsing read strings. Do not use
   * without external parser.
   * @param fl
   * File with aligned reads fasta format
   * @return
   * Iterable collection of reads
   */
  def initFastaReads(fl: File): Iterable[Read] = {
    val seqs = readFastaDNASequence(fl)
    val readsMap = flip(seqs.map(en => (en._1, en._2.getSequenceAsString)).toMap)
    val samRecords = toSAMRecords(readsMap)
    val reads = toReads(samRecords)
    initReadFreqs(reads, readsMap)
    log("Number of distinct reads: %d".format(reads.size))
    reads
  }

  /**
   * Read alignment in internal format.
   * Reads in lines aligned with spaces
   * all have the same length. Output
   * of alignment postprocessing tool.
   * (Temporary solution)
   * @param fl
   * Aligned reads (Internal txt format)
   * @return
   * Iterable collection of reds
   */
  @deprecated
  def initTxtReads(fl: File): Iterable[Read] = {
    val lines = fromFile(fl).getLines
    val readsMap = flip(lines.zipWithIndex.map(s => ("Read" + s._2, s._1)).toMap)
    val samRecords = toSAMRecords(readsMap)
    val reads = toReads(samRecords)
    initReadFreqs(reads, readsMap.map(entry => (entry._1, entry._2)).toMap)
    reads
  }

  /**
   * Init read frequencies according to counter map
   * @param reads
   * Collection of reads
   * @param readsMap
   * Counter Map
   */
  private def initReadFreqs(reads: Iterable[Read], readsMap: Map[String, Set[String]]) {
    reads foreach (r => {
      r.freq = readsMap(r.seq).size;
      r.ids = readsMap(r.seq)
    })
  }

  /**
   * Cover SAMRecords into Read objects
   * @param sams
   * Collection of SAMRecords
   * @return
   * Collection of reads
   */
  private def toReads(sams: Iterable[SAMRecord]) = {
    for (sam <- sams) yield new Read(sam)
  }

  /**
   * Cover Map Entries with string and SAMRecord
   * into Read objects
   * @param sams
   * Map with strings and SAMRecords
   * @return
   * Collection of Reads
   */
  private def toReads(sams: Map[String, SAMRecord]) = {
    for (e <- sams) yield {
      e._2.setReadString(e._1)
      e._2.setAlignmentStart(1)
      new Read(e._2)
    }
  }

  /**
   * Converts list of Strings into counter map:
   * i. e. (A,A,B,C,C,C) -> ({A:2},{B,1},{C:3})
   * @param lines
   * Iterable collection of strings
   * @return
   * Map with counts of strings
   */
  private def toCounterMap(lines: Iterator[String]): mutable.Map[String, Int] = {
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

  /**
   * Flip map into mulimap
   * Map[X, Y] => Map[ Y, Set[X] ]
   * @param m
   * Input Map[X, Y]
   * @tparam X
   * Type of Key
   * @tparam Y
   * Type of Value
   * @return
   * Map[ Y, Set[X] ]
   */
  def flip[X, Y](m: Map[X, Y]): Map[Y, Set[X]] =
    m.groupBy(_._2).mapValues(_.map(_._1).toSet)

  def log(mes: String) = {
    val date = new Date()
    println("%s %s".format(sdf.format(date), mes))
  }

  /**
   * Deserialize SAMRecords from strings
   * @param reads
   * Reads in strings
   * @return
   * SAMRecords collection
   */
  @deprecated
  private def toSAMRecords(reads: Map[String, Set[String]]) = {
    for (r <- reads)
    yield toSAMRecord(r)
  }

  /**
   * Deserialize one {@see SAMRecord} from {@see String}
   * @param st
   * { @see SAMRecord} in { @see String}
   * @return
   * Wrapped { @see SAMRecord} object
   */
  @deprecated
  private def toSAMRecord(st: (String, Set[String])) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(1)
    sam.setReadString(st._1)
    sam
  }

  protected[exec] def printGreetings = {
    log(LINE)
    log(KGEM_STR.format(Main.getClass.getPackage.getImplementationVersion))
    log(LINE)
  }
}
