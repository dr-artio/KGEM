package edu.gsu.cs.kgem

import collection.mutable
import io.SAMParser
import model.Read
import util.Random
import java.io.File
import scala.io.Source._
import net.sf.samtools.{SAMFileHeader, SAMRecord}

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/30/13
 * Time: 3:02 PM
 * To change this template use File | Settings | File Templates.
 */
package object exec {

  /**
   * Read SAM file and fold according to extended
   * sequences.
   * @param fl
   *           SAM file
   * @return
   *         Iterable collection of Reads
   */
  def initSAMReads(fl: File): Iterable[Read] = {
    val samRecords = SAMParser.readSAMFile(fl)
    var extSAMRecords = samRecords.map(s => SAMParser.toExtendedString(s))
    val l = extSAMRecords.map(s => s.length).max
    extSAMRecords = samRecords.map(s => SAMParser.toExtendedString(s, l))
    val readsMap = toCounterMap(extSAMRecords.iterator)
    val samMap = samRecords.map(s => (SAMParser.toExtendedString(s, l), s)).toMap
    val reads = toReads(samMap)
    initReadFreqs(reads, readsMap)
    reads
  }

  /**
   * Old method for parsing read strings. Do not use
   * without external parser.
   * @param fl
   *           File wit extended reads in extended format
   * @return
   *         Iterable collection of reads
   */
  @deprecated
  def initTXTReads(fl: File): Iterable[Read] = {
    val lines = fromFile(fl).getLines
    val readsMap = toCounterMap(lines)
    val samRecords = toSAMRecords(readsMap.keys)
    val reads = toReads(samRecords)
    initReadFreqs(reads, readsMap)
    reads
  }

  /**
   * Init read frequencies according to counter map
   * @param reads
   *              Collection of reads
   * @param readsMap
   *                 Counter Map
   */
  private def initReadFreqs(reads: Iterable[Read], readsMap: mutable.Map[String, Int]) {
    reads foreach (r => {
      r.freq = readsMap(r.seq)
    })
  }

  /**
   * Cover SAMRecords into Read objects
   * @param sams
   *             Collection of SAMRecords
   * @return
   *         Collection of reads
   */
  private def toReads(sams: Iterable[SAMRecord]) = {
    for (sam <- sams) yield new Read(sam)
  }

  /**
   * Cover Map Entries with string and SAMRecord
   * into Read objects
   * @param sams
   *             Map with strings and SAMRecords
   * @return
   *         Collection of Reads
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
   *              Iterable collection of strings
   * @return
   *         Map with counts of strings
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
   * Deserialize SAMRecords from strings
   * @param reads
   *              Reads in strings
   * @return
   *         SAMRecords collection
   */
  @deprecated
  private def toSAMRecords(reads: Iterable[String]) = {
    for (r <- reads)
    yield toSAMRecord(r)
  }

  /**
   * Deserialize one SAMRecord from String
   * @param st
   * @return
   */
  @deprecated
  private def toSAMRecord(st: String) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(1)
    sam.setReadString(st)
    sam
  }
}
