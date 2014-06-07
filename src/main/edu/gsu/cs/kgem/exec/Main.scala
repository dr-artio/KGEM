package edu.gsu.cs.kgem.exec


/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/15/13
 * Time: 4:02 PM
 * Main entry point of the program.
 * Arguments passed through command line
 */
object Main {
  /**
   * Main method for console
   * run of kGEM tool.
   * System.exit() called at the end.
   * Please use doMain if called as
   * sub procedure
   * @param args
   *             Arguments
   */
  def main(args: Array[String]): Unit = {
    doMain(args)
    sys.exit(0)
  }

  /**
   * Method for execution main procedure.
   * Does not call System.exit and good
   * for calls from java
   * @param args
   *             Arguments
   */
  def doMain(args: Array[String]): Unit = {
    val s = System.currentTimeMillis
    parseArgs(args)
    printGreetings()
    log("Reading input files...")
    initInputData()
    log("KGEM started...")
    val gens = executeKgem()
    log("Output results...")
    outputResults(gens, s)
  }
}
