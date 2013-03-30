package edu.gsu.cs.kgem.io

import net.sourceforge.argparse4j.ArgumentParsers
import java.io.{OutputStream, FileOutputStream, PrintStream, File}
import net.sourceforge.argparse4j.inf.ArgumentParserException

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 2:58 PM
 * To change this template use File | Settings | File Templates.
 */
object ArgumentParser {

  private val READS_PARAMETER = "reads"
  private val K_PARAMETER = "k"
  private val OUTPUT_PARAMETER = "out"

  /**
   * Method for parsing command line parameters
   * and handle mistakes in parameters list
   * @param args
   *             Array of command line parameters
   * @return
   *         Tuple of parsed and converted parameters
   */
  def parse(args: Array[String]) = {
    val parser = ArgumentParsers.newArgumentParser("KGEM")
      .description("Error correction based on KGEM.")

    var k = 50
    var fl: File = null
    var out = System.out

    parser.addArgument(READS_PARAMETER)
      .metavar("reads file")
      .help("File containing preprocessed sequencing data"
      + " file with extension (.txt) or (.sam) "
      + "reads in extended format")
      .`type`(classOf[File])

    parser.addArgument("-k").dest(K_PARAMETER)
      .metavar("sample size")
      .`type`(classOf[Integer])
      .help("Parameter k - the size of sample being randomly chosen "
      + "as seeds. Depends on expectation of variability expected "
      + "on exploring region (Default: " + k + ")")

    parser.addArgument("-o", "-out").dest(OUTPUT_PARAMETER)
      .metavar("output filename")
      .setDefault[PrintStream](out)
      .`type`(classOf[FileOutputStream])
      .help("Output file name. (Default: stdout)")

    try {
      val n = parser.parseArgs(args)
      val kk = n.get(K_PARAMETER).asInstanceOf[Int]
      k = if (kk > 1) kk else k
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
