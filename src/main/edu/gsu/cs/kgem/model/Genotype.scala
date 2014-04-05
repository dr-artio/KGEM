package edu.gsu.cs.kgem.model

import org.biojava3.core.sequence.compound.{DNACompoundSet, NucleotideCompound}
import collection.JavaConversions._
import collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/16/13
 * Time: 7:26 PM
 * To change this template use File | Settings | File Templates.
 */
object Genotype {
  /**
   * Initialization of mappings for building one-to-one
   * correspondence between NucleotideCompounds and
   * string representations of nucleotides
   * (current DNA alphabet)
   */
  val nuclMap: Map[NucleotideCompound, String] = {
    DNACompoundSet.getDNACompoundSet.getAllCompounds.map(n => {
      var s = n.getShortName.toUpperCase
      if (s == "N") s = "-"
      (n, s)
    }).toMap
  }
  val sMap: Map[String, Set[NucleotideCompound]] = reverseMap(nuclMap)

  val N = "N"
  val eps = 0.0025
  private var id = 0

  private def generateID = {
    id += 1
    id
  }

  def reverseMap[K, V](mp: Map[K, V]): Map[V, Set[K]] = {
    val r = mp.values.map(v => (v, mutable.Set[K]())).toMap
    mp.foreach((e: (K, V)) => {
      r(e._2) += e._1
    })
    return r.map((e: (V, mutable.Set[K])) => e._1 -> e._2.toSet)
  }

  def nameForCompound(n: NucleotideCompound) = nuclMap(n)
}

import Genotype._

class Genotype(n: Int) {
  val ID = generateID

  val data: List[mutable.Map[String, Double]] = {
    List.range(0, n).map(i => mutable.Map(sMap.keys.map(s => (s, 0D)).toSeq: _*))
  }
  var freq = 1.0
  var convergen = false

  def this(str: String) = {
    this(str.length)
    var i = 0
    for (d <- data) {
      val cur_symb = str(i).toString
      if (d.contains(cur_symb))
        d(cur_symb) = 1.0
      else {
        val avg = 1.0 / d.size
        d.keys.foreach(k => d(k) = avg)
      }
      i += 1
    }
    round
  }

  def this(reads: List[Read]) = {
    this(reads.map(r => r.end).max)
    for (r <- reads) addRead(r)
    round
  }

  /**
   * Add read to the genotype, i. e. increase
   * the the value on each position covered by
   * read and corresponding nucleotide
   * @param r
   * Read object with alignment information
   * and sequence.
   */
  @inline
  def addRead(r: Read) = {
    var s = r.beg
    r.seq foreach (c => {
      data(s)(c.toString) += 1
      s += 1
    })
  }

  def sqNormalize = {
    normalize
    data foreach (m => {
      val s = m.values.map(v => v * v).sum
      if (s != 0) m foreach (e => m(e._1) /= s)
    })
  }

  private def normalize = {
    data foreach (m => {
      val s = m.values.sum
      if (s != 0) m foreach (e => m(e._1) /= s)
    })
  }

  def round = {
    data foreach (m => {
      val ss = (m.map(x => x._2).sum * 1.1) / m.size
      val s = m.maxBy(e => e._2)
      m foreach (e => m(e._1) = eps)
      if (s._2 > ss) m(s._1) = 1
    })
    sqNormalize
    data
  }

  def toIntegralString = {
    val s = new StringBuilder
    data foreach (m => s ++= {
      val avg = 1.01 * m.map(e => e._2).sum / m.size
      val mm = m.maxBy(e => e._2)
      if (mm._2 > avg) mm._1
      else N
    })
    s.toString
  }

  override def toString = {
    val s = new StringBuilder
    data foreach (m => s ++= (m + "\n"))
    s.toString
  }
}
