package edu.gsu.cs.kgem.model

import org.biojava3.core.sequence.DNASequence

/**
 * Created by aartsiomenka1 on 3/28/2014.
 */
object Clusterer {

  /**
   * Cluster original reads
   * @param gens
   * Collection of Haplotypes
   * @param reads
   * Collection of reads
   * @param pqrs
   * Clustering coefficients
   */
  def getClusteredFastaSequencies(gens: Iterable[Genotype], reads: Iterable[Read],
                                  pqrs: Array[Array[Double]], rs: Map[String, DNASequence], k: Int) = {
    var i = -1
    val genSeqs = gens.map(g => (g, g.toIntegralString.replaceAll("-", "").replaceAll("N", ""))).toList

    val groups = cluster(reads, pqrs, genSeqs, rs)

    groups.map(rd => {
      i += 1
      rd._2.map(r => {
        val ds = rs(r._1.ids.head)
        ds.setOriginalHeader("c%d_h%d_%s_%d".format(rd._1.ID, rd._1.ID, ds.getOriginalHeader, r._1.ids.size))
        ds
      })
    }).flatten
  }

  private def cluster(reads: Iterable[Read], pqrs: Array[Array[Double]],
                      genSeqs: List[(Genotype, String)],
                      rs: Map[String, DNASequence]) = {
    val res = finalizeClustering(reads, pqrs, Map[(Read, Int), Genotype](), genSeqs, rs)

    res
  }

  private def finalizeClustering(reads: Iterable[Read], pqrs: Array[Array[Double]],
                                 gbReads: Map[(Read, Int), Genotype],
                                 genSeqs: List[(Genotype, String)],
                                 rs: Map[String, DNASequence]) = {
    reads.zipWithIndex.filter(x => !gbReads.contains(x) || gbReads(x) != null).groupBy(rd => {
      if (gbReads.contains(rd)) {
        gbReads(rd)
      } else {
        val cl = column(pqrs, rd._2)
        val index = cl.indexWhere(p => p == cl.max)
        genSeqs(index)._1
      }
    }).toMap
  }

  private def column[A, M[_]](matrix: M[M[A]], colIdx: Int)
                             (implicit v1: M[M[A]] => Seq[M[A]], v2: M[A] => Seq[A]): Seq[A] =
    matrix.map(_(colIdx))


  /**
   * Merge two maps with collections as
   * values. Result map will contain union
   * of keys as new keyset and joint lists of
   * values
   * @param m1
   * Operand 1
   * @param m2
   * Operand 2
   * @tparam K
   * Generic parameter type of Key
   * @tparam V
   * Generic parameter type of Value
   * @return
   */
  private def merge[K, V](m1: Map[K, Iterable[V]], m2: Map[K, Iterable[V]]): Map[K, Iterable[V]] = {
    val k1 = Set(m1.keysIterator.toList: _*)
    val k2 = Set(m2.keysIterator.toList: _*)
    val intersection = k1 & k2

    val r1 = for (key <- intersection) yield (key -> (m1(key) ++ m2(key)))
    val r2 = m1.filterKeys(!intersection.contains(_)) ++ m2.filterKeys(!intersection.contains(_))
    r2 ++ r1
  }

  private def mean[T <% Double](vals: Iterable[T]) = {
    val N = vals.size
    vals.map(_.toDouble).sum / N
  }

  private def sigma[T <% Double](vals: Iterable[T]) = {
    val meanValue = mean(vals)
    Math.sqrt(mean(vals.map(x => Math.pow(x - meanValue, 2))))
  }
}
