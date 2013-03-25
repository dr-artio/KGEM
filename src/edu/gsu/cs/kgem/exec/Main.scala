package edu.gsu.cs.kgem.exec

import collection.mutable
import edu.gsu.cs.kgem.model.{KGEM, Genotype, Read}
import net.sf.samtools.{SAMFileHeader, SAMRecord, SAMFileReader}
import scala.io.Source._
import scala.util.Random

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/15/13
 * Time: 4:02 PM
 * To change this template use File | Settings | File Templates.
 */
object Main {

  def main(args: Array[String]) {
    val s = System.currentTimeMillis
    val lines = fromFile("grinder-reads.txt").getLines
    val readsMap = toCounterMap(lines)
    val SAMRecords = toSAMRecords(readsMap.keys)
    val reads = toReads(SAMRecords)
    initReadFreqs(reads, readsMap)
    val seeds = sample(reads, 50)
    var gens = seeds.map(s => new Genotype(s.seq)).toList
    KGEM.initReads(reads.toList)
    gens = KGEM.run(gens)
    for (g <- gens) {
      println(g.freq + "\n" + g.toIntegralString.replace("-", ""))
    }
    println("The whole procedure took " +
      ((System.currentTimeMillis - s) * 0.0001 / 6) + " minutes\n" +
      "Total number of haplotypes is " + gens.size + "\nbye bye")
  }

  def sample[T](iter: Iterable[T], size: Int) = {
    if (iter.size < size) iter.toList
    var res = new mutable.MutableList[T]()
    val rnd = new Random(System.currentTimeMillis)
    var needed = size
    var len = iter.size
    val iterator = iter.iterator
    while (needed > 0 && iterator.hasNext) {
      val item = iterator.next
      if (rnd.nextInt(len) < needed) {
        res += item
        needed -= 1
      }
      len -= 1
    }
    res

  }

  def initReadFreqs(reads: Iterable[Read], readsMap: mutable.Map[String, Int]) {
    reads foreach (r => {
      r.freq = readsMap(r.seq)
    })
  }

  def toReads(sams: Iterable[SAMRecord]) = {
    for (sam <- sams) yield new Read(sam)
  }

  def toCounterMap(lines: Iterator[String]): mutable.Map[String, Int] = {
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

  def toSAMRecords(reads: Iterable[String]) = {
    for (r <- reads)
    yield toSAMRecord(r)
  }

  def toSAMRecord(st: String) = {
    val sam = new SAMRecord(new SAMFileHeader)
    sam.setAlignmentStart(0)
    sam.setReadString(st)
    sam
  }
}
