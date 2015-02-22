package edu.gsu.cs.kgem.model

import org.apache.commons.math3.distribution.BinomialDistribution
import org.biojava3.core.sequence.compound.{DNACompoundSet, NucleotideCompound}
import collection.JavaConversions._
import collection.mutable
import scala.collection.mutable.ListBuffer

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
  val N = 'N'
  var eps = 0.0025
  private var id = 0


  val nuclMap: Map[NucleotideCompound, String] = {
    DNACompoundSet.getDNACompoundSet.getAllCompounds.map(n => {
      var s = n.getShortName.toUpperCase
      if (s.charAt(0) == N) s = "-"
      (n, s)
    }).toMap
  }
  val sMap: Map[String, Set[NucleotideCompound]] = reverseMap(nuclMap)

  private def generateID = {
    id += 1
    id
  }

  def reverseMap[K, V](mp: Map[K, V]): Map[V, Set[K]] = {
    val r = mp.values.map(v => (v, mutable.Set[K]())).toMap
    mp.foreach((e: (K, V)) => {
      r(e._2) += e._1
    })
    r.map((e: (V, mutable.Set[K])) => e._1 -> e._2.toSet)
  }

  def nameForCompound(n: NucleotideCompound) = nuclMap(n)
}

import Genotype._

class Genotype(n: Int) {
  private var p_value = 0.05
  val ID = generateID
  private var epsi = -1.0
  private val erronius_snps = new ListBuffer[(Char, Int)]()
  var size = 0.0
  val coverage = new Array[Double](n)
  val reads = new ListBuffer[Read]()
  val table = Array.fill[mutable.Map[Char, ListBuffer[Read]]](n) {
     mutable.Map(sMap.keys.map(s => s.charAt(0) -> new ListBuffer[Read]()).toSeq: _*)
  }

  val data: List[mutable.Map[Char, Double]] = {
    List.range(0, n).map(i => mutable.Map(sMap.keys.map(s => (s.charAt(0), 0D)).toSeq: _*))
  }
  var freq = 1.0
  var convergen = false

  def this(str: String) = {
    this(str.length)
    addRead(str)
    round
  }

  def this(reads: Iterable[Read]) = {
    this(reads.map(r => r.end).max)
    for (r <- reads) addRead(r)
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
  def addRead(r: Read): Unit = {
    addRead(r.seq, r.beg, r.freq)
    reads += r
    addReadToTable(r)
  }

  @inline
  def removeRead(r: Read): Unit = {
    if (reads contains r) {
      removeRead(r.seq, r.beg, r.freq)
      reads -= r
      removeReadFromTable(r)
    }
  }

  private def addReadToTable(r: Read): Unit = {
    for (c <- r.seq.zipWithIndex.par) {
      if (table(c._2) contains c._1)
        table(c._2)(c._1) += r
    }
  }

  private def removeReadFromTable(r: Read): Unit = {
    for (c <- r.seq.zipWithIndex.par) {
      if (table(c._2) contains c._1)
        table(c._2)(c._1) -= r
    }
  }

  private def addRead(str: String, b: Int = 0, freq: Double = 1.0): Unit = {
    size += freq
    epsi = -1.0
    str.zipWithIndex.par foreach (c => {
      val s = c._2
      val d = data(s)
      val cur_symb = c._1
      if (d.contains(cur_symb)) {
        coverage(s) += freq
        d(cur_symb) += freq
      }
    })
  }

  private def removeRead(str: String, b: Int = 0, freq: Double = 1.0): Unit = {
    size -= freq
    epsi = -1.0
    str.zipWithIndex.par foreach (c => {
      val s = c._2
      val d = data(s)
      val cur_symb = c._1
      if (d.contains(cur_symb)) {
        if (coverage(s) > 0.0) {
          coverage(s) -= freq
          d(cur_symb) -= freq
        }
      }
    })
  }

  def getSecondHaplotype = {
    val snps = findCorrelatedSNPs
    if (snps.size > 1) {
      val majorHaplotype = toIntegralString
      val secondHaplotype = new StringBuffer()
      for (c <- majorHaplotype.zipWithIndex) {
        val snp = snps.find(s => s._2 == c._2)
        if (snp != None)
          secondHaplotype.append(snp.get._1)
        else
          secondHaplotype.append(c._1)
      }
      secondHaplotype.toString
    } else ""
  }

  def getSecondGenotype = {
    val snps = findCorrelatedSNPs
    if (snps.size > 1) {
      val sreads = reads.filter(r => snps.forall(s => r.seq(s._2) == s._1))
      if (sreads.size > 2) {
        //p_value += 0.05
        new Genotype(sreads)
      } else
        None
    }  else {
      None
    }
  }

  def findCorrelatedSNPs: Iterable[(Char, Int)] = {
    val snps = data.zipWithIndex
//      .filter(m => {
//      val s = coverage(m._2)
//      if (s > 0)
//        (1 - m._1.values.max / s) > epsilon
//      else
//        false
//    })
    var selectedSnps = snps.map(s => {
      val M = s._1.maxBy(_._2)
      val m = s._1.filter(_ != M).maxBy(_._2)
      (m._1, s._2, m._2)
    }).sortBy(-_._3).takeWhile(x => x._3 >= KGEM.threshold)

    val correlatedSnps = new ListBuffer[(Char, Int)]()

    while(correlatedSnps.size < 2 && selectedSnps.size > 1) {
      correlatedSnps.clear()
      val firstSnp = selectedSnps.head
      correlatedSnps += ((firstSnp._1, firstSnp._2))
      val snps = selectedSnps.par.filter(x => linkageDisequilibrium(x._2, x._1, firstSnp._2, firstSnp._1))

      if (snps.nonEmpty){
        val snp = snps.maxBy(x => correlation(x._2, x._1, firstSnp._2, firstSnp._1)._2)
        correlatedSnps += ((snp._1, snp._2))
        return correlatedSnps.toList
      } else {
        selectedSnps = selectedSnps.tail
      }
    }

    //correlatedSnps.clear()
    correlatedSnps.toList
  }

  def linkageDisequilibrium(i: Int, ci: Char, j: Int, cj: Char): Boolean = {
    val x = correlation(i, ci, j, cj)
    x._1 <= p_value / Math.pow(table.length, 2) * 2
  }

  def intersect(A: Iterable[Read], B: Iterable[Read]): Iterable[Read] = {
    val AB = A ++ B
    val sAB = AB.toList.sortBy(_.seq)

    sAB.sliding(2).filter(x => x.head == x.tail.head).map(_.head).toList
  }

  def correlation(i: Int, ci: Char, j: Int, cj: Char): (Double, Double) = {

    // Default return
    val d = (1.0, .0)

    if (Math.abs(i-j) <= 1 || (ci == '-' && cj == '-')) return d
    val MA = data(i).maxBy(_._2)
    val MB = data(j).maxBy(_._2)
    val x_11 = intersect(table(i)(MA._1), table(j)(MB._1))
    val x_12 = intersect(table(i)(MA._1), table(j)(cj))
    val x_21 = intersect(table(i)(ci), table(j)(MB._1))
    val x_22 = intersect(table(i)(ci), table(j)(cj))
    var X_11 = x_11.map(_.freq).sum
    var X_12 = x_12.map(_.freq).sum
    var X_21 = x_21.map(_.freq).sum
    var X_22 = x_22.map(_.freq).sum
    val X = X_11 + X_12 + X_21 + X_22
    val p = Math.max((X_12 * X_21 / X_11) / X, 10 / X)
    if (p > 1) return d
    val dist = new BinomialDistribution(X.toInt, p)
//    X_11 /= X
//    X_12 /= X
//    X_21 /= X
//    X_22 /= X

    val res = 1.0 - dist.cumulativeProbability(X_22.toInt)
    (res, X_22)
  }


  def normalizedValue(i: Int, c: Char): Double = {
    if (data.length > i && i >= 0 && data(i).contains(c))
      return data(i)(c) / coverage(i)
    0
  }

  def sqNormalize() = {
    normalize()
    data foreach (m => {
      val s = m.values.map(v => v * v).sum
      if (s != 0) m foreach (e => m(e._1) /= s)
    })
  }

  private def normalize() = {
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
    sqNormalize()
    data
  }

  def toIntegralString = {
    val s = new StringBuilder
    data foreach (m => s.append({
      val avg = 1.001 * m.map(_._2).sum / m.size
      val mm = m.maxBy(_._2)
      if (mm._2 > avg) mm._1
      else N
    }))
    s.toString()
  }

  def epsilon = {
    if (epsi < 0) {
      val s = data.map(m => {
        val s = m.values.sum
        if (s > 0) 1 - m.values.max / s else 0
      }).sum
      epsi = s / data.length
    }
    epsi
  }

  override def toString = {
    val s = new StringBuilder
    data foreach (m => s ++= (m + "\n"))
    s.toString()
  }
}
