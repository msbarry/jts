/*
 * Copyright (c) 2021 Martin Davis.
 *
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License 2.0
 * and Eclipse Distribution License v. 1.0 which accompanies this distribution.
 * The Eclipse Public License is available at http://www.eclipse.org/legal/epl-v20.html
 * and the Eclipse Distribution License is available at
 *
 * http://www.eclipse.org/org/documents/edl-v10.php.
 */
package org.locationtech.jts.index;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.index.chain.MonotoneChain;
import org.locationtech.jts.math.MathUtil;
import org.locationtech.jts.util.IntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * A semi-static spatial index for points which occur 
 * in a spatially-coherent sequence.
 * In particular, this is suitable for indexing the vertices
 * of a {@link LineString} or {@link Polygon} ring.
 * <p>
 * The index is constructed in a batch fashion on a given sequence of coordinates.
 * <p>
 * Note that this index queries only the individual points
 * of the input coordinate sequence, 
 * <b>not</b> any line segments which might be lie between them.
 * <p>
 * The input coordinate array is read-only, 
 * and is not changed when vertices are removed.
 * 
 * @author Martin Davis
 *
 */
public class MonotoneChainSequencePackedRtree {

  /**
   * Number of items/nodes in a parent node.
   * Determined empirically.  Performance is not too sensitive to this.
   */
  private static final int NODE_CAPACITY = 16;
  public static final MonotoneChain[] EMPTY = new MonotoneChain[0];
  private final double overlap;

  private MonotoneChain[] items;
  private double[] nodeBounds;
  private double[] itemBounds;
  private int[] levelOffset;
  private int nodeCapacity  = NODE_CAPACITY;

  /**
   * Creates a new tree over the given sequence of coordinates.
   * The sequence should be spatially coherent to provide query performance.
   *
   * @param mcs a sequence of points
   */
  public MonotoneChainSequencePackedRtree(List<MonotoneChain> mcs, double overlap) {
    this.items = mcs.toArray(EMPTY);
    this.overlap = overlap;
    build();
  }
  
  private void build() {
    levelOffset = computeLevelOffsets();
    itemBounds = new double[items.length * 4];
    int pos = 0;
    for (MonotoneChain item : items) {
      Envelope env = item.getEnvelope(overlap);
      itemBounds[pos++] = env.getMinX();
      itemBounds[pos++] = env.getMinY();
      itemBounds[pos++] = env.getMaxX();
      itemBounds[pos++] = env.getMaxY();
    }
    nodeBounds = createBounds();
  }

  /**
   * Computes the level offsets.
   * This is the position in the <tt>bounds</tt> array of each level.
   * 
   * The levelOffsets array includes a sentinel value of offset[0] = 0.
   * The top level is always of size 1,
   * and so also indicates the total number of bounds.
   * 
   * @return the level offsets
   */
  private int[] computeLevelOffsets() {
    IntArrayList offsets = new IntArrayList();
    offsets.add(0);
    int levelSize = items.length;
    int currOffset = 0;
    do {
      levelSize = levelNodeCount(levelSize);
      currOffset += levelSize;
      offsets.add(currOffset);
    } while (levelSize > 1);
    return offsets.toArray();
  }

  private int levelNodeCount(int numNodes) {
    return MathUtil.ceil(numNodes, nodeCapacity);
  }
  
  private double[] createBounds() {
    int boundsSize = levelOffset[levelOffset.length - 1] + 1;
    double[] bounds = new double[boundsSize * 4];
    Arrays.fill(bounds, Double.NaN);
    fillItemBounds(bounds);
    
    for (int lvl = 1; lvl < levelOffset.length; lvl++) {
      fillLevelBounds(lvl, bounds);
    }
    return bounds;
  }
  
  private void fillLevelBounds(int lvl, double[] bounds) {
    int levelStart = levelOffset[lvl - 1]; 
    int levelEnd = levelOffset[lvl];
    int nodeStart = levelStart;
    int levelBoundIndex = levelOffset[lvl];
    do {
      int nodeEnd = MathUtil.clampMax(nodeStart + nodeCapacity, levelEnd);
      computeEnvelope(bounds, bounds, levelBoundIndex++, nodeStart, nodeEnd);
      nodeStart = nodeEnd;
    }
    while (nodeStart < levelEnd);
  }

  private void fillItemBounds(double[] bounds) {
    int nodeStart = 0;
    int boundIndex = 0;
    do {
      int nodeEnd = MathUtil.clampMax(nodeStart + nodeCapacity, items.length);
      computeEnvelope(itemBounds, bounds, boundIndex++, nodeStart, nodeEnd);
      nodeStart = nodeEnd;
    }
    while (nodeStart < items.length);
  }
  
  private static void computeEnvelope(double[] in, double[] out, int outIndex, int start, int end) {
    double minX = Double.POSITIVE_INFINITY;
    double minY = Double.POSITIVE_INFINITY;
    double maxX = Double.NEGATIVE_INFINITY;
    double maxY = Double.NEGATIVE_INFINITY;
    int pos = start * 4;
    for (int i = start; i < end; i++) {
      minX = Math.min(minX, in[pos++]);
      minY = Math.min(minY, in[pos++]);
      maxX = Math.max(maxX, in[pos++]);
      maxY = Math.max(maxY, in[pos++]);
    }
    out[outIndex * 4] = minX;
    out[outIndex * 4 + 1] = minY;
    out[outIndex * 4 + 2] = maxX;
    out[outIndex * 4 + 3] = maxY;
  }
  
  //------------------------

  /**
   * Queries the index to find all items which intersect an extent.
   * The query result is a list of the indices of input coordinates
   * which intersect the extent.
   * 
   * @param queryEnv the query extent
   * @return an array of the indices of the input coordinates
   */
  public List<MonotoneChain> query(Envelope queryEnv) {
    List<MonotoneChain> resultList = new ArrayList<>();
    int level = levelOffset.length - 1;
    queryNode(queryEnv, level, 0, resultList);
    return resultList;
  }


  private static boolean intersects(double[] bounds, int nodeIndex, Envelope env) {
    boolean isBeyond = (env.getMaxX() < bounds[nodeIndex])
            || (env.getMaxY() < bounds[nodeIndex+1])
            || (env.getMinX() > bounds[nodeIndex+2])
            || (env.getMinY() > bounds[nodeIndex+3]);
    return ! isBeyond;
  }
  
  private void queryNode(Envelope queryEnv, int level, int nodeIndex, List<MonotoneChain> resultList) {
    int boundsIndex = (levelOffset[level] + nodeIndex) * 4;
    //--- node is empty
    if (Double.isNaN(nodeBounds[boundsIndex]))
      return;
    if (intersects(nodeBounds, boundsIndex, queryEnv)) {
      int childNodeIndex = nodeIndex * nodeCapacity;
      if (level == 0) {
        queryItemRange(queryEnv, childNodeIndex, resultList);
      } else {
        queryNodeRange(queryEnv, level - 1, childNodeIndex, resultList);
      }
    }
  }

  private void queryNodeRange(Envelope queryEnv, int level, int nodeStartIndex, List<MonotoneChain> resultList) {
    int levelMax = levelSize(level);
    for (int i = 0; i < nodeCapacity; i++) {
      int index = nodeStartIndex + i;
      if (index >= levelMax) 
        return;
      queryNode(queryEnv, level, index, resultList);
    }    
  }

  private int levelSize(int level) {
    return levelOffset[level + 1] - levelOffset[level];
  }

  private void queryItemRange(Envelope queryEnv, int itemIndex, List<MonotoneChain> resultList) {
    for (int i = 0; i < nodeCapacity; i++) {
      int index = itemIndex + i;
      if (index >= items.length) 
        return;
      if (intersects(itemBounds, index * 4, queryEnv))
        resultList.add(items[index]);
    }
  }
}
