package edu.gsu.cs.kgem.io

import org.scalatest.FlatSpec
/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/18/13
 * Time: 12:18 PM
 */
/**
 *  User: nicholas
 *  Date: 7/18/13
 *  Time: 12:18 PM
 */
class ArgumentParserTest extends FlatSpec {

  "The parseArguments function" should " return empty Option[Config]" +
    " when no arguments are passed." in {
    val args = new Array[String](0)
    val config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)
  }

  it should " return empty Option[Config] when bad k is supplied" in {
    val args = new Array[String](2)
    args(0) = "reads.sam"
    args(1) = "0"
    var config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "-5"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

  }

  it should " return empty Option[Config] when bad threshold is supplied"  in {
    val args = new Array[String](3)
    args(0) = "reads.sam"
    args(1) = "-t"
    args(2) = "0"
    var config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "-t"
    args(2) = "-5"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)
  }

  it should "return empty Option[Config] when bad range is supplied" in {
    val args = new Array[String](2)
    args(0) = "reads.sam"
    args(1) = "-10:15"
    var config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "0:15"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "15:0"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "0:0"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "-10:15"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "0:15"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "15:0"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "0:0"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "10:15:20"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "funfun"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)
  }

  it should " return non-empty Option[Config] when a good k is supplied" in {
    val args = new Array[String](2)
    args(0) = "reads.sam"
    args(1) = "10"
    var config = ArgumentParser.parseArguments(args)
    assert(!config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "10:15"
    config = ArgumentParser.parseArguments(args)
    assert(!config.isEmpty)
  }

  it should "return empty Option[Config] when a bad scoring function is specified" in {
    val args = new Array[String](3)
    args(0) = "reads.sam"
    args(1) = "-f"
    args(2) = "TIC"
    var config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)

    args(0) = "reads.sam"
    args(1) = "--scoring-func"
    args(2) = "TIC"
    config = ArgumentParser.parseArguments(args)
    assert(config.isEmpty)
  }
}
