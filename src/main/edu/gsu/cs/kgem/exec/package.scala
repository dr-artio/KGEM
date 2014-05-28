package edu.gsu.cs.kgem

import edu.gsu.cs.kgem.model._
import collection.mutable
import edu.gsu.cs.kgem.io.{ArgumentParser, SAMParser}
import java.io.{PrintStream, File}
import net.sf.samtools.{SAMFileHeader, SAMRecord}
import org.biojava3.core.sequence.io.FastaReaderHelper.readFastaDNASequence
import collection.JavaConversions._
import scala.io.Source.fromFile
import java.util.Date
import java.text.SimpleDateFormat
import edu.gsu.cs.kgem.io.OutputHandler._
import org.biojava3.core.sequence.DNASequence
import edu.gsu.cs.kgem.io.Config
import scala.Some
import scala.collection.parallel.{ForkJoinTasks, ForkJoinTaskSupport, TaskSupport}
import scala.concurrent.forkjoin.ForkJoinPool

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:02 PM
 */
package object exec {
  val USER_DIR = "user.dir"
  val KGEM_STR = "kGEM version %s: Local Reconstruction for Mixed Viral Populations."
  var numproc: TaskSupport = null

  //private values and variables
  private val sdf = new SimpleDateFormat("[hh:mm:ss a]")
  private val FASTA = Array[String](".fas", ".fa", ".fasta")
  private val LINE = "-----------------------------------------------------------"
  private var out: PrintStream = null
  private var config = new Config()
  private var seqs: List[DNASequence] = null
  private var reads: List[Read] = null
  private var k: Int = -1
  private var threshold: Int = 0
  private var n: Int = 0
  private var seeds: Iterable[Genotype] = null

  //Public methods
  /**
   * Method to perform KGEM. When call
   * from another program as library
   * all parameters are mandatory
   * @param reads
   * Collection od reads { @see Read}
   * @param k
   * Parameter k in model (number of initial candidates)
   * @param threshold
   * distance threshold for initial guesses
   * @param eps
   * Error rate (will be taken quarter of value actually sent)
   * should be below 0.5 and greater than 0
   * @param pr_threshold
   * Frequency threshold for dropping rare sequences
   * @param seeds
   * Initial seeds or null if not given
   * @return
   * Set of genotypes corresponding to a given set of reads
   */
  def executeKgem(reads: List[DNASequence] = seqs, k: Int = k, numproc: TaskSupport = config.numproc, threshold: Int = threshold,
                  eps: Double = config.epsilon, pr_threshold: Double = config.prThr,
                  seeds: Iterable[Genotype] = seeds): List[Genotype] = {
    if (numproc != null ) {
      this.numproc = numproc
      log("Numprocs set %d.".format(this.numproc.parallelismLevel))
      setParallelismGlobally(this.numproc.parallelismLevel)
    }
    this.reads = convertFastaReads(reads).toList
    n = this.reads.map(r => r.freq).sum.toInt
    log("Pre")
    KGEM.initReads(this.reads.toList)
    log("Pre")
    Genotype.eps = eps
    val gens = if (seeds == null) {

      if (pr_threshold >= 0) KGEM.initThreshold(pr_threshold)
      else KGEM.initThreshold()

      val seeds = MaxDistanceSeedFinder.findSeeds(this.reads, k, threshold)
      KGEM.run(seeds)
    } else {
      KGEM.run(seeds)
    }
    gens.view.toIndexedSeq.sortBy(-_.freq).toList
  }

  // Protected section
  protected[exec] def parseArgs(args: Array[String]) = {
    ArgumentParser.parseArguments(args) match {
      case None => sys.exit(1)
      case Some(config: Config) =>
        this.config = config
        setupOutput(config.output) match {
          case None => sys.exit(1)
          case Some((out: PrintStream)) => this.out = out
        }
    }
  }

  protected[exec] def setParallelismGlobally(numThreads: Int): Unit = {
    val parPkgObj = scala.collection.parallel.`package`
    val defaultTaskSupportField = parPkgObj.getClass.getDeclaredFields.find{
      _.getName == "defaultTaskSupport"
    }.get

    defaultTaskSupportField.setAccessible(true)
    defaultTaskSupportField.set(
      parPkgObj,
      new scala.collection.parallel.ForkJoinTaskSupport(
        new scala.concurrent.forkjoin.ForkJoinPool(numThreads)
      )
    )
  }

  protected[exec] def initInputData(readsFile: File = config.readsFile) = {
    if (!FASTA.exists(readsFile.getName.toLowerCase.endsWith)) sys.exit(1)
    seqs = readFastaDNASequence(readsFile).values.toList
    k = config.k
    threshold = config.threshold
    if (config.consensusFile != null)
      seeds = initFastaReads(config.consensusFile).map(r => new Genotype(r.seq))
  }

  protected[exec] def outputResults(gens: List[Genotype], s: Long) = {
    def str_modifier: (String => String) = if (config.is_cleaned)
      s => s.replaceAll("-", "")
    else
      s => s
    if (config.is_reads)
      outputResult(out, gens, n, str_modifier)
    else
      outputHaplotypes(out, gens, str_modifier)
    log("The whole procedure took %.2f minutes".format((System.currentTimeMillis - s) * 0.0001 / 6))
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
  protected[exec] def initSAMReads(fl: File): Iterable[Read] = {
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
  protected[exec] def initFastaReads(fl: File): Iterable[Read] = {
    val seqs = readFastaDNASequence(fl).values.toList
    convertFastaReads(seqs)
  }

  /**
   * Transform DNASequence reads to
   * internal kgem objects
   * @param seqs
   * Collection of DNASequence objects
   * @return
   * Collection of Read objects
   */
  private def convertFastaReads(seqs: Iterable[DNASequence]): Iterable[Read] = {
    val readsMap = flip(seqs.map(en => (en.getOriginalHeader, en.getSequenceAsString)).toMap)
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
  protected[exec] def initTxtReads(fl: File): Iterable[Read] = {
    val lines = fromFile(fl).getLines()
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
      r.freq = readsMap(r.seq).size
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
  private def toSAMRecords(reads: Map[String, Set[String]]) = {
    for (r <- reads)
    yield toSAMRecord(r)
  }

  /**
   * Deserialize one SAMRecord from String
   * @param st
   * { @see SAMRecord} in { @see String}
   * @return
   * Wrapped { @see SAMRecord} object
   */
  private def toSAMRecord(st: (String, Set[String])) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(1)
    sam.setReadString(st._1)
    sam
  }

  protected[exec] def printGreetings() = {
    log(LINE)
    log(KGEM_STR.format(Main.getClass.getPackage.getImplementationVersion))
    log(LINE)
  }
}
