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
  def main(args: Array[String]): Unit = {
    val s = System.currentTimeMillis
    parseArgs(args)
    printGreetings
    log("Reading input files...")
    initInputData()
    log("KGEM started...")
    val gens = executeKgem()
    log("Output results...")
    outputResults(gens, s)
  }
}
