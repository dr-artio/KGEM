package edu.gsu.cs.kgem.io

import java.io.File
import scopt.OptionParser
import edu.gsu.cs.kgem.exec.{KGEM_STR, Main}
import scala.collection.parallel.{ForkJoinTaskSupport, TaskSupport}
import scala.concurrent.forkjoin.ForkJoinPool

case class Config(readsFile: File = null, k: Int = 50, threshold: Int = 3,
                  consensusFile: File = null, prThr: Double = -1, epsilon: Double = 0.0025,
                  is_reads: Boolean = false, is_cleaned: Boolean = false,
                  output: File = null, numproc: TaskSupport = null)


object ArgumentParser {

  def parseArguments(args: Array[String]): Option[Config] = {
    val parser = new OptionParser[Config]("kGEM") {
      head(KGEM_STR.format(Main.getClass.getPackage.getImplementationVersion))
      arg[File]("ReadsFile") action {
        (x, c) => c.copy(readsFile = x)
      } text "Fasta (*.fa, *.fas, *.fasta) file containing aligned sequence data."
      arg[Int]("k") action {
        (x, c) => c.copy(k = x)
      } validate {
        case x => if (x > 0) success
        else failure("Must be positive" +
          " integer value.")
      } text "Number of initial seeds for clustering."
      opt[Double]('t', "frequency threshold") action {
        (x, c) => c.copy(prThr = x)
      } validate {
        x => if (x >= 0 && x <= 0.5) success else failure("Frequency threshold should be 0<=x<=0.5")
      } text "Frequency threshold"
      opt[Double]('e', "epsilon") action {
        (x, c) => c.copy(epsilon = x)
      } validate {
        x => if (x >= 0 && x <= 0.25) success else failure("Error rate should be 0<=x<=0.25")
      } text "Approximate flat error rate"
      opt[Int]('d', "threshold") action {
        (x, c) => c.copy(threshold = x)
      } validate {
        x => if (x > 0) success else failure("Threshold must be greater than 0.")
      } text "Min Hamming distance between seeds for init stage."
      opt[Int]('n', "numproc") action {
        (x, c) => c.copy(numproc = new ForkJoinTaskSupport(new ForkJoinPool(x)))
      } text "Number of threads for parallel execution"
      opt[File]('g', "seeds-file") action {
        (x, c) => c.copy(consensusFile = x)
      } text "Optional Fasta File containing initial seeds."
      opt[File]('o', "output-filename") action {
        (x, c) => c.copy(output = x)
      } validate {
        x => if (x == null) failure("Output file must be specified. Use -o <filepath>")
        else success
      } text "File to output results."
      opt[Unit]('r', "out-as-reads") action {
        (x, c) => c.copy(is_reads = true)
      } text "Use if reads as separate file required (frequencies through counts)"
      opt[Unit]('u', "unaligned") action {
        (x, c) => c.copy(is_cleaned = true)
      } text "Output unaligned (removed alignment information)"
      opt[Unit]('h', "help") action {
        (x, c) => showUsage; sys.exit(0)
      } text "Prints this help text."
      opt[Unit]('v', "version") action {
        (x, c) => showHeader; sys.exit(0)
      } text "Prints kGEM version."
    }
    // parser.parse returns Option[C]
    parser.parse(args, Config())
  }
}
