package edu.gsu.cs.kgem.test

import org.scalatest.FunSuite
import scala.collection.mutable
import edu.gsu.cs.kgem.model.MaxDistanceWrapper.{hammingDistance, min}

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 4/13/13
 * Time: 3:20 PM
 * To change this template use File | Settings | File Templates.
 */
class TestSuite extends FunSuite {

  test("test Hamming distance method") {
    val str1 = "ACCG-TC"
    val str2 = "ACT-TTC"
    val result = hammingDistance(str1, str2)
    assert(result == 1)
  }
  test("Test min method") {
    val i = 1
    val j = 11
    assert(min(i, j) == i)
    assert(min(j, i) == i)
  }

  test("Complicated test for population") {
    val pop = new mutable.MutableList[String]()
    pop += "A-----ACTG"
    pop += "-A----ACCC"
    pop += "--A---AAAA"
    pop += "---A--CCTG"
    pop += "----A-ACTT"
    pop += "-----ACCTC"
    val k = 6
    val t = 2
    val first = pop(0)
    var seeds = new mutable.MutableList[String]()
    seeds += first
    var distanceMap = pop.filter(r => !r.equals(first)).map(r => (r, hammingDistance(first, r)))
    while (seeds.size < k) {
      val cur = distanceMap.maxBy(e => e._2)
      //if (cur._2 < t) break
      println("Current max HD: %d".format(cur._2))
      println(cur)
      seeds += cur._1
      distanceMap = distanceMap.map(e => (e._1, min(e._2, hammingDistance(cur._1, e._1)))) //.filter(e => (e._2 >= threshold))
      println("Current size of DistnceMap: %d".format(distanceMap.size))
      //distanceMap = readArr.filter(r => !seeds.contains(r)).map(r => (r, seeds.map(s => hammingDistance(s, r)).min)).toMap
    }
    println(seeds)
  }
}
