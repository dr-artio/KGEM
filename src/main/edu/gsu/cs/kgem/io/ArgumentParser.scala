package edu.gsu.cs.kgem.io

import java.io.File
import scopt.OptionParser
import edu.gsu.cs.kgem.exec._

case class Config(readsFile: File = null, k: Int = 50, threshold: Int = 3, clustering: File = null,
                  consensusFile: File = null, prThr: Double = 0, output: File = new File(System.getProperty(USER_DIR)))


object ArgumentParser {

  private val _singlePattern = """(\A\d+\Z)""".r
  private val _doublePattern = """(\A\d+:\d+\Z)""".r

  def parseArguments(args: Array[String]): Option[Config] = {
    val parser = new OptionParser[Config]("kGEM") {
      head(KGEM_STR.format(Main.getClass.getPackage.getImplementationVersion))
      arg[File]("ReadsFile") action {
        (x, c) => c.copy(readsFile = x)
      } text ("Fasta (*.fa, *.fas, *.fasta) file containing aligned sequence data.")
      arg[Int]("k") action {
        (x, c) => c.copy(k = x)
      } validate {
        case x => if (x > 0) success
        else failure("Must be positive" +
          " integer value.")
      } text ("Number of initial seeds for clustering.")
      opt[Double]('t', "frequency threshold") action {
        (x, c) => c.copy(prThr = x)
      } text ("Frequency threshold")
      opt[Int]('d', "threshold") action {
        (x, c) => c.copy(threshold = x)
      } validate {
        x => if (x > 0) success else failure("Threshold must be greater than 0.")
      } text ("Min Hamming distance between seeds for init stage.")
      opt[File]('g', "consensus-file") action {
        (x, c) => c.copy(consensusFile = x)
      } text ("Optional Fasta File containing initial seeds.")
      opt[File]('o', "output-dir") action {
        (x, c) => c.copy(output = x)
      } text ("Directory to output results.")
      opt[File]('c',"clustering") action { (x, c) => c.copy(clustering = x)}
      opt[Unit]('h', "help") action {
        (x, c) => showUsage; sys.exit(0)
      } text ("Prints this help text.")
      opt[Unit]('v', "version") action {
        (x, c) => showHeader; sys.exit(0)
      } text ("Prints kGEM version.")
    }
    // parser.parse returns Option[C]
    return parser.parse(args, Config())
  }

  private def isValidRange(arg: String): Boolean = arg match {
    case this._doublePattern(c) => {
      val rng = arg.split(":").map(x => x.toInt)
      return rng(0) > 0 && rng(1) >= rng(0)
    }
    case this._singlePattern(c) => {
      return arg.toInt > 0
    }
    case _ => return false
  }

  private def parseRange(arg: String): Range.Inclusive = arg match {
    case this._doublePattern(c) => {
      val rng = arg.split(":").map(x => x.toInt)
      (rng(0) to rng(1))
    }
    case this._singlePattern(c) => {
      return (arg.toInt to arg.toInt)
    }
    case _ => return (0 to 0)
  }
}
