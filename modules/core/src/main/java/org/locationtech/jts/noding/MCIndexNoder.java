/*
 * Copyright (c) 2016 Vivid Solutions, and others.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
package org.locationtech.jts.noding;

import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateXY;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.index.SegmentSequencePackedRtree;
import org.locationtech.jts.index.SpatialIndex;
import org.locationtech.jts.index.chain.MonotoneChain;
import org.locationtech.jts.index.chain.MonotoneChainOverlapAction;
import org.locationtech.jts.index.hprtree.HPRtree;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * Nodes a set of {@link SegmentString}s using a index based
 * on {@link MonotoneChain}s and a {@link SpatialIndex}.
 * The {@link SpatialIndex} used should be something that supports
 * envelope (range) queries efficiently (such as a <code>Quadtree</code>}
 * or {@link HPRtree} (which is the default index provided).
 * <p>
 * The noder supports using an overlap tolerance distance .
 * This allows determining segment intersection using a buffer for uses
 * involving snapping with a distance tolerance.
 *
 * @version 1.7
 */
public class MCIndexNoder
    extends SinglePassNoder
{
  private final List<SegmentString> segmentStrings = new ArrayList<>();
  private final SpatialIndex<SegmentSequencePackedRtree> index = new HPRtree<>();
  private int idCounter = 0;
  private Collection nodedSegStrings;
  private double overlapTolerance = 0;

  public MCIndexNoder()
  {
  }
  
  public MCIndexNoder(SegmentIntersector si)
  {
    super(si);
  }

  /**
   * Creates a new noder with a given {@link SegmentIntersector}
   * and an overlap tolerance distance to expand intersection tests with.
   * 
   * @param si the segment intersector
   * @param overlapTolerance the expansion distance for overlap tests
   */
  public MCIndexNoder(SegmentIntersector si, double overlapTolerance)
  {
    super(si);
    this.overlapTolerance = overlapTolerance;
  }

  public List getMonotoneChains() { return segmentStrings; }

  public SpatialIndex getIndex() { return index; }

  public Collection getNodedSubstrings()
  {
    return  NodedSegmentString.getNodedSubstrings(nodedSegStrings);
  }

  public void computeNodes(Collection inputSegStrings)
  {
//    System.err.println("computeNodes(" + inputSegStrings + ")");
    this.nodedSegStrings = inputSegStrings;
    for (Object inputSegString : inputSegStrings) {
      SegmentString segStr = (SegmentString) inputSegString;
      Envelope env = new Envelope();
      for (Coordinate coord : segStr.getCoordinates()) {
        env.expandToInclude(coord);
      }
      env.expandBy(overlapTolerance);
      index.insert(env, new SegmentSequencePackedRtree(segStr, segmentStrings.size()));
      segmentStrings.add(segStr);
    }
    Envelope envelope = new Envelope();
    int id = 0;
    for (SegmentString queryChain : segmentStrings) {
//      System.err.println("  queryChain=" + queryChain);
      for (int node1 = 0; node1 < queryChain.size() - 1; node1++) {
//        System.err.println("    node1=" + node1 + " " + new CoordinateXY(queryChain.getCoordinate(node1)) + " " + new CoordinateXY(queryChain.getCoordinate(node1 + 1)));
        envelope.setToNull();
        envelope.expandToInclude(queryChain.getCoordinate(node1));
        envelope.expandToInclude(queryChain.getCoordinate(node1 + 1));
        envelope.expandBy(overlapTolerance);
        for (SegmentSequencePackedRtree item : index.query(envelope)) {
//          System.err.println("      testChain=" + item.getSegmentString());
          if (item.getId() >= id) {
            for (int node2 : item.query(envelope)) {
//              System.err.println("        node2=" + node2 + " " + new CoordinateXY(item.getSegmentString().getCoordinate(node2)) + " " + new CoordinateXY(item.getSegmentString().getCoordinate(node2 + 1)));
              segInt.processIntersections(queryChain, node1, item.getSegmentString(), node2);
              // short-circuit if possible
              if (segInt.isDone())
                return;
            }
          }
        }
      }
      id++;
    }

//    for (int i = 0; i < segmentStrings.size(); i++) {
//      SegmentString queryChain = segmentStrings.get(i);
//      for (int j = i; j < segmentStrings.size(); j++) {
//        SegmentString testChain = segmentStrings.get(j);
//        for (int node1 = 0; node1 < queryChain.size() - 1; node1++) {
//          for (int node2 = 0; node2 < testChain.size() - 1; node2++) {
//            segInt.processIntersections(queryChain, node1, testChain, node2);
//            // short-circuit if possible
//            if (segInt.isDone())
//              return;
//          }
//        }
//      }
//    }
//    intersectChains();
//System.out.println("MCIndexNoder: # chain overlaps = " + nOverlaps);
  }

  private void intersectChains()
  {

//    Envelope envelope = new Envelope();
//    int id = 0;
//    for (SegmentString queryChain : monoChains) {
//      for (int node1 = 0; node1 < queryChain.size() - 1; node1++) {
//        envelope.setToNull();
//        envelope.expandToInclude(queryChain.getCoordinate(node1));
//        envelope.expandToInclude(queryChain.getCoordinate(node1 + 1));
//        for (SegmentSequencePackedRtree item : index.query(envelope)) {
//          if (item.getId() > id) {
//            for (int node2 : item.query(envelope)) {
//              segInt.processIntersections(queryChain, node1, item.getSegmentString(), node2);
//              nOverlaps++;
//              // short-circuit if possible
//              if (segInt.isDone())
//                return;
//            }
//          }
//        }
//      }
//      id++;
//    }


    for (int i = 0; i < segmentStrings.size(); i++) {
      SegmentString queryChain = segmentStrings.get(i);
      for (int j = i + 1; j < segmentStrings.size(); j++) {
        SegmentString testChain = segmentStrings.get(j);
        for (int node1 = 0; node1 < queryChain.size() - 1; node1++) {
          for (int node2 = 0; node2 < testChain.size() - 1; node2++) {
            segInt.processIntersections(queryChain, node1, testChain, node2);
            // short-circuit if possible
            if (segInt.isDone())
              return;
          }
        }
      }
    }
  }

//  private void add(SegmentString segStr)
//  {
//    monoChains.add(segStr);
//    Envelope env = new Envelope();
//    for (Coordinate coord : segStr.getCoordinates()) {
//      env.expandToInclude(coord);
//    }
//    env.expandBy(overlapTolerance);
//    index.insert(env, new SegmentSequencePackedRtree(segStr, idCounter++));
//  }
}
