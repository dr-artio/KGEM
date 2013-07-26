package edu.gsu.cs.kgem.io

import java.io.File
import net.sf.samtools.{CigarElement, SAMFileReader, SAMRecord}
import collection.mutable
import collection.JavaConversions._

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/28/13
 * Time: 7:29 PM
 * To change this template use File | Settings | File Templates.
 */
object SAMParser {
  val S = " "

  def readSAMFile(f: File): Iterable[SAMRecord] = {
    if (!f.exists()) System.err.println("File not found!")
    val reader = new SAMFileReader(f)
    var res = new mutable.MutableList[SAMRecord]()
    val iter = reader.iterator
    while (iter.hasNext) res += iter.next
    res
  }

  def toExtendedString(sam: SAMRecord, l: Int = 0): String = {
    val sb = new StringBuilder
    var index = 0
    val str = sam.getReadString
    val cs = sam.getCigar.getCigarElements
    for (i <- 0 until sam.getAlignmentStart - 1)
      sb.append(S)
    for (cigar: CigarElement <- cs) {
      if (!cigar.getOperator.consumesReadBases) {
        for (i <- 0 until cigar.getLength)
          sb.append("-")
      }
      else {
        for (i <- 0 until cigar.getLength) {
          sb.append(str(index))
          index += 1
        }
      }
    }
    val cutStr = sb.toString
    var ll = cutStr.length
    if (ll >= l)
      return cutStr
    else
      while (ll < l) {
        ll += 1
        sb.append(S)
      }
    sb.toString
  }

}
