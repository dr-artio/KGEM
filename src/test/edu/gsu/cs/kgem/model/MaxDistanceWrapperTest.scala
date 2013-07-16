package edu.gsu.cs.kgem.model

import org.scalatest.FlatSpec
import edu.gsu.cs.kgem.model.MaxDistanceWrapper.{hammingDistance, min}

/**
 *  Author: Nicholas Mancuso
 */
class MaxDistanceWrapperTest extends FlatSpec {
  "The Hamming distance function" should "sum the number of non-space or non-dash " +
    "mismatches between two equal-length strings" in {
    val str1 = "ACCG-TC"
    val str2 = "ACT-TTC"
    val result = hammingDistance(str1, str2)
    assert(result == 1)
  }

  it should "throw IllegalArgumentException if the arguments are of unequal length" in  {
    val str1 = "ACCG"
    val str2 = "ACT-TTC"
    intercept[IllegalArgumentException] {
      hammingDistance(str1, str2)
    }
  }

  "The Min function" should "return the minimum between two integers" in {
    val i = 1
    val j = 11
    assert(min(i, j) == i)
    assert(min(j, i) == i)
  }
}