package edu.gsu.cs.kgem.model.selection

import org.scalatest.FlatSpec
import edu.gsu.cs.kgem.model.selection.Score._

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/17/13
 * Time: 4:56 PM
 * To change this template use File | Settings | File Templates.
 */
class ScoreTest extends FlatSpec {
  "The getScoringFunction function" should "return the requested function" +
    "when the name is valid." in {
    var name = Score.AIC_
    var func = getScoringFunction(name)
    assert(func != null)

    name = Score.AICc_
    func = getScoringFunction(name)
    assert(func != null)

    name = Score.BIC_
    func = getScoringFunction(name)
    assert(func != null)

    name = Score.BICMAP_
    func = getScoringFunction(name)
    assert(func != null)
  }

  it should "throw an IllegalArgumentException when the name is invalid" in {
    val name = "FOO"
    intercept[IllegalArgumentException] {
      getScoringFunction(name)
    }
  }
}
