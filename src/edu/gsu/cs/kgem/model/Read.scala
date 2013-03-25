package edu.gsu.cs.kgem.model

import net.sf.samtools.SAMRecord

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/17/13
 * Time: 10:55 PM
 * To change this template use File | Settings | File Templates.
 */
class Read(rc: SAMRecord) {
  var freq = 1.0
  val beg = rc.getAlignmentStart
  val seq = rc.getReadString
  val len = rc.getReadLength
  val end = beg + len

  override def equals(obj: Any) = {
    if (obj.isInstanceOf[Read]) {
      val or = obj.asInstanceOf[Read]
      seq.equals(or.seq) && beg == or.beg
    }
    false
  }
}
