package edu.gsu.cs.kgem.io

import java.io.File
import scopt.OptionParser

case class Config(readsFile: File = null, k: Range.Inclusive = (50 to 50), threshold: Int = 3,
  scoringFunc: String = "AICc", consensusFile: File = null, prThr: Double = -1, clustering: File = null,
  output: File = new File("./"))


object ArgumentParser {

  private val _singlePattern = """(\A\d+\Z)""".r
  private val _doublePattern = """(\A\d+:\d+\Z)""".r

  def parseArguments(args: Array[String]): Option[Config] = {
    val parser = new OptionParser[Config]("kGEM") {
      head("kGEM version 0.4.2: Local Reconstruction for Mixed Viral Populations.")
      arg[File]("ReadsFile") action {
        (x, c) => c.copy(readsFile = x)
      } text("Fasta or Sam file containing aligned sequence data.")
      arg[String]("k") action {
        (x, c) => c.copy(k = parseRange(x))
      } validate {
        case x => if (isValidRange(x)) success else failure("Lower bound must" +
          " be greater than 0 and upper bound must be greater than or equal to lower" +
          " bound.")
      } text("Number of initial seeds for clustering. May be a single positive\n" +
        "\tinteger value or may be a range to try for model selection (e.g. 5:10).")
      opt[Double]('t', "frequency threshold") action {
        (x, c) => c.copy(prThr = x)
      } text("Frequency threshold")
      opt[Int]('d', "threshold") action {
        (x, c) => c.copy(threshold = x)
      } validate {
        x => if (x > 0) success else failure("Threshold must be greater than 0.")
      } text("Min Hamming distance between seeds for init stage.")
      opt[String]('f', "scoring-func") action {
        (x, c) => c.copy(scoringFunc = x)
      } validate {
        x => if (List("AIC", "AICc", "BIC", "AICc", "BICMAP").contains(x)) success
             else failure("Scoring function must be AIC, AICc, BIC, or BICMAP.")
      } text("The scoring method to use in model selection. Only used when a\n" +
        "\trange is specified. (Default: \'AICc\').\n" +
        "\tAccepted Values:\n\n" +
        "\t\tAIC: Akaike Information Criteria\n" +
        "\t\tAICc: Corrected Akaike Information Criteria\n" +
        "\t\tBIC: Bayesian Information Criteria\n" +
        "\t\tBICMAP: Bayesian Information Criteria using MAP estimate")
      opt[File]('g', "consensus-file") action {
        (x, c) => c.copy(consensusFile = x)
      } text("Optional Fasta File containing initial seeds.")
      opt[File]('o', "output-dir") action {
        (x, c) => c.copy(output = x)
      } text("Directory to output results.")
      opt[File]('c',"clustering") action { (x, c) => c.copy(clustering = x)}
      opt[Unit]('h', "help") action { (x, c) => showUsage; sys.exit(0) } text("Prints this help text.")
      opt[Unit]('v', "version") action { (x, c) => showHeader; sys.exit(0) } text("Prints kGEM version.")
      }
    // parser.parse returns Option[C]
    return parser.parse(args, Config())
  }

  private def isValidRange(arg:String): Boolean = arg match {
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
