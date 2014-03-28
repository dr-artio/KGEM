package edu.gsu.cs.kgem.model

import net.sf.samtools.SAMRecord

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/17/13
 * Time: 10:55 PM
 * Wrapper for Read
 */
class Read(rc: SAMRecord) {
  var freq = 1.0
  val beg = rc.getAlignmentStart - 1
  val seq = rc.getReadString
  val len = rc.getReadLength
  val end = beg + len
  var ids: Set[String] = null

  override def equals(obj: Any) = {
    if (obj.isInstanceOf[Read]) {
      val or = obj.asInstanceOf[Read]
      seq.equals(or.seq) && beg == or.beg
    }
    false
  }
}
