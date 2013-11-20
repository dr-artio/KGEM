package edu.gsu.cs.kgem

import collection.mutable
import io.SAMParser
import model.Read
import java.io.File
import net.sf.samtools.{SAMFileHeader, SAMRecord}
import org.biojava3.core.sequence.io.FastaReaderHelper.readFastaDNASequence
import collection.JavaConversions._
import scala.io.Source.fromFile
import org.biojava3.core.sequence.DNASequence

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
   * SAM file
   * @return
   * Iterable collection of Reads
   */
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
    println("Number of distinct reads: %d".format(reads.size))
    reads
  }

  /**
   * Read alignment in internal format.
   * Reads in lines aligned with spaces
   * all have the same length. Output
   * of alignment postprocessing tool.
   * (Temporary solution)
   * @param fl
   *           Aligned reads (Internal txt format)
   * @return
   *         Iterable collection of reds
   */
  def initTxtReads(fl: File): Iterable[Read] = {
    val lines = fromFile(fl).getLines
    val readsMap = flip(lines.zipWithIndex.map(s => ("Read"+s._2, s._1)).toMap)
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
   *          Input Map[X, Y]
   * @tparam X
   *           Type of Key
   * @tparam Y
   *           Type of Value
   * @return
   *         Map[ Y, Set[X] ]
   */
  def flip[X, Y](m: Map[X, Y]): Map[Y, Set[X]] =
    m.groupBy(_._2).mapValues(_.map(_._1).toSet)

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
   *           {@see SAMRecord} in {@see String}
   * @return
   *         Wrapped {@see SAMRecord} object
   */
  @deprecated
  private def toSAMRecord(st: (String, Set[String])) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(1)
    sam.setReadString(st._1)
    sam
  }
}
