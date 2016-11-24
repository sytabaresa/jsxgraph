/*
    Copyright 2008-2016
        Matthias Ehmann,
        Darko Draculic,
        Michael Gerhaeuser,
        Carsten Miller,
        Alfred Wassermann

    This file is part of JSXGraph.

    JSXGraph is free software dual licensed under the GNU LGPL or MIT License.

    You can redistribute it and/or modify it under the terms of the

      * GNU Lesser General Public License as published by
        the Free Software Foundation, either version 3 of the License, or
        (at your option) any later version
      OR
      * MIT License: https://github.com/jsxgraph/jsxgraph/blob/master/LICENSE.MIT

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License and
    the MIT License along with JSXGraph. If not, see <http://www.gnu.org/licenses/>
    and <http://opensource.org/licenses/MIT/>.
 */

 /*global JXG: true, define: true*/
 /*jslint nomen: true, plusplus: true*/

 /* depends:
  jxg
  math/math
  utils/type
  */

 /**
  * @fileoverview In this file the conic sections defined.
  */

define([
     'jxg', 'math/math', 'math/geometry', 'utils/type'
], function (JXG, Mat, Geometry, Type) {

    "use strict";

    Mat.Implicitplot = {
        is_contour_present: function(f, x, y, d) {
            var fxy = f(x, y),
                fxdy = f(x + d, y),
                fxyd = f(x, y + d),
                fxdyd = f(x + d, y + d);

            if (fxy * fxdy < 0 ||
                fxy * fxyd < 0 ||
                fxdy * fxdyd < 0) {
                return true;
            } else if (fxy === 0 &&
                        (fxdy * fxdyd < 0 ||
                         fxyd * fxdyd < 0 ||
                         fxdyd === 0)) {
                return true;
            } else if (fxdy === 0 &&
                        (fxy * fxyd < 0 ||
                         fxyd * fxdyd < 0 ||
                         fxyd === 0)) {
                return true;
            }
            return false;
        },

        subdivide_rectangle: function(f, array, depth, x, y, d, searchDepth) {
            this.find_points(f, array, depth + 1, x,           y,           d * 0.5, searchDepth);
            this.find_points(f, array, depth + 1, x + d * 0.5, y,           d * 0.5, searchDepth);
            this.find_points(f, array, depth + 1, x + d * 0.5, y + d * 0.5, d * 0.5, searchDepth);
            this.find_points(f, array, depth + 1, x,           y + d * 0.5, d * 0.5, searchDepth);
        },

        find_points: function(f, array, depth, x, y, d, searchDepth) {
            var firstDepth = searchDepth,
                secondDepth = searchDepth + 2;

            if (depth < firstDepth) {
                this.subdivide_rectangle(f, array, depth, x, y, d, searchDepth);
            } else if (this.is_contour_present(f, x, y, d)) {
                if (depth < secondDepth) {
                    this.subdivide_rectangle(f, array, depth, x, y, d, searchDepth);
                } else {
                    this.addToArray(array,
                            [(2 * x + d) * 0.5,
                             (2 * y + d) * 0.5]);
                }
            }
        },

        addToArray: function(ptsArray, point) {
            var pos, len;

            if (ptsArray.length === 0) {
                ptsArray[0] = [point];
            } else {
                len = ptsArray.length;
                for (pos = 0; pos < len && ptsArray[pos][0][0] < point[0]; pos++) {}

                if (pos == len) {
                    ptsArray[pos] = [point];
                } else if (point[0] != ptsArray[pos][0][0]) {
                    ptsArray.splice(pos, 0, [point]);
                }
                //ptsArray[pos].push(point);
            }
        },

        function organizePoints(f, ptsArray, step) {
            var sortedPtsArray,
                curvePart,
                stepX, stepY,
                lastPoint, indexOfLastPoint,
                noPoints, i, j, len,
                iter,
                neighbours, previousIndex, nextIndex;
			var pos;
			var noPossiblePoints, minYin, maxYin = 0, startIndex = 0;
			var neighbourIndex, lastIndex;
			var noX, noY, sX, sY, ptsSX, signSX;
			var next, result, parts;
			var minDistance, idx, isReverse, back;
            var d;

            if (ptsArray.length === 0) {
                return [];
            }

            sortedPtsArray = [];
            curvePart = 0;
            // find curve parts

            stepX = step;
            stepY = step;

            sortedPtsArray[curvePart] = [ptsArray[0][0]];
            lastPoint = ptsArray[0][0];
            indexOfLastPoint = 0;

            ptsArray[0].splice(0, 1);
            if (ptsArray[0].length == 0) {
                ptsArray.splice(0, 1);
            }

            noPoints = 0;
            len = ptsArray.length;
            for (i = 0; i < len; i++) {
                for (j = 0; j < ptsArray[i].length; j++) {
                    noPoints++;
                }
            }

            for (iter = 0; iter < noPoints; iter++) {
                neighbours = [];

                previousIndex = (indexOfLastPoint != 0 &&
                        ptsArray[indexOfLastPoint - 1].length != 0 &&
                        equal(ptsArray[indexOfLastPoint-1][0][0] + stepX, lastPoint[0])) ?
                            indexOfLastPoint - 1 : -1;
                nextIndex = (indexOfLastPoint != ptsArray.length - 1 &&
                        ptsArray[indexOfLastPoint+1].length != 0 &&
                        equal(ptsArray[indexOfLastPoint+1][0][0] - stepX, lastPoint[0]) ) ?
                            indexOfLastPoint + 1 : -1;

                for (i = 0; i < ptsArray[indexOfLastPoint].length; i++) {
                    if (equal(lastPoint[1], ptsArray[indexOfLastPoint][i][1]+stepY) ||
                        equal(lastPoint[1], ptsArray[indexOfLastPoint][i][1]-stepY)) {
                        neighbours.push([indexOfLastPoint, i]);
                    }
                }

                if (previousIndex != -1) {
                    for (i = 0; i < ptsArray[previousIndex].length; i++) {
                        if (equal(lastPoint[1], ptsArray[previousIndex][i][1])) {
                            neighbours.push([previousIndex, i]);
                        }
                    }
                }

                if (nextIndex != -1) {
                    for(i = 0; i < ptsArray[nextIndex].length; i++) {
                        if (equal(lastPoint[1], ptsArray[nextIndex][i][1])) {
                            neighbours.push([nextIndex, i]);
                        }
                    }
                }

                if (neighbours.length === 0) {
                    for (pos = 0; pos < ptsArray.length; pos++) {
                        if (ptsArray[pos].length != 0) {
                            break;
                        }
                    }
                    if (pos == ptsArray.length) {
                        break;
                    }

                    noPossiblePoints = ptsArray[pos].length;
                    minYin = 0;
                    maxYin = 0;
                    startIndex = 0;

                    if (noPossiblePoints > 1) {
                        neighbourIndex = 0;
                        for (i = 1; i < noPossiblePoints; i++) {
                            if (ptsArray[pos][i][1] < ptsArray[pos][minYin][1]) {
                                minYin = i;
                                maxYin = i;
                            }
                        }

                        for (i = 0; i < noPossiblePoints; i++) {
                            if (equal(ptsArray[pos][maxYin][1], ptsArray[pos][i][1] - stepY)) {
                                maxYin = i;
                            }
                        }

                        if (pos + 1 < ptsArray.length - 1 &&
                            ptsArray[pos + 1].length > 1) {
                            for (i = 0; i < ptsArray[pos + 1].length; i++) {
                                if (equal(ptsArray[pos][minYin][1], ptsArray[pos+1][i][1])) {
                                    neighbourIndex = minYin;
                                    break;
                                }
                                if (equal(ptsArray[pos][maxYin][1], ptsArray[pos+1][i][1])) {
                                    neighbourIndex = maxYin;
                                    break;
                                }
                            }
                        }
                        startIndex = (neighbourIndex == maxYin) ? minYin : maxYin;
                    }
                    sortedPtsArray[++curvePart] = [];
                    sortedPtsArray[curvePart].push(ptsArray[pos][startIndex]);

                    lastPoint = ptsArray[pos][startIndex];
                    indexOfLastPoint = pos;
                    ptsArray[pos].splice(startIndex,1);
                } else {
                    if (sortedPtsArray[curvePart].length == 1) { // start of closed curve
                        sortedPtsArray[curvePart].push(ptsArray[neighbours[0][0]][neighbours[0][1]]);
                        lastPoint = ptsArray[neighbours[0][0]][neighbours[0][1]];
                        indexOfLastPoint = neighbours[0][0];
                        ptsArray[neighbours[0][0]].splice(neighbours[0][1],1);
                    } else if (neighbours.length == 1) { // just one neighbour point
                        sortedPtsArray[curvePart].push(ptsArray[neighbours[0][0]][neighbours[0][1]]);
                        lastPoint = ptsArray[neighbours[0][0]][neighbours[0][1]];
                        indexOfLastPoint = neighbours[0][0];
                        ptsArray[neighbours[0][0]].splice(neighbours[0][1],1);

                        if (sortedPtsArray[curvePart].length > 3) {
                            lastIndex = sortedPtsArray[curvePart].length-1;
                            if ( (  equal(sortedPtsArray[curvePart][lastIndex-2][0], sortedPtsArray[curvePart][lastIndex-1][0]) &&
                                    equal(sortedPtsArray[curvePart][lastIndex-1][0], sortedPtsArray[curvePart][lastIndex][0] )) ||
                                 (  equal(sortedPtsArray[curvePart][lastIndex-2][1], sortedPtsArray[curvePart][lastIndex-1][1]) &&
                                    equal(sortedPtsArray[curvePart][lastIndex-1][1], sortedPtsArray[curvePart][lastIndex][1]) ) ) {
                                        sortedPtsArray[curvePart].splice(lastIndex-1, 1);
                            }
                        }
                    } else if (neighbours.length > 1) {
                        noX = 0;
                        noY = 0;
                        ptsSX = [];

                        if (equal(lastPoint[1], ptsArray[neighbours[0][0]][neighbours[0][1]][1])) {
                            sX = ptsArray[neighbours[0][0]][neighbours[0][1]][0] - lastPoint[0];
                        } else {
                            sX = ptsArray[neighbours[1][0]][neighbours[1][1]][0] - lastPoint[0];
                        }
                        if (equal(lastPoint[0], ptsArray[neighbours[0][0]][neighbours[0][1]][0])) {
                            sY = ptsArray[neighbours[0][0]][neighbours[0][1]][1] - lastPoint[1];
                        } else {
                            sY = ptsArray[neighbours[1][0]][neighbours[1][1]][1] - lastPoint[1];
                        }

                        ptsSX[0] = [];
                        ptsSX[1] = [];

                        signSX = (sX > 0) ? 1 : -1;

                        for(i = indexOfLastPoint + signSX; i < ptsArray.length && i >= 0; i += signSX) {
                            noX = 0;
                            for(j = 0; j < ptsArray[i].length; j++) {
                                if (equal(lastPoint[1], ptsArray[i][j][1])) {
                                    noX++;
                                    ptsSX[0].push([i, j]);
                                }
                                if (equal(lastPoint[1]+sY, ptsArray[i][j][1])) {
                                    ptsSX[1].push([i, j]);
                                }
                            }
                            if (noX == 0) {
                                break;
                            }
                        }

                        for (i = 0; i < ptsSX[0].length / 2; i++) {
                            sortedPtsArray[curvePart].push(ptsArray[ptsSX[0][i][0]][ptsSX[0][i][1]]);
                            lastPoint = ptsArray[ptsSX[0][i][0]][ptsSX[0][i][1]];
                            indexOfLastPoint = ptsSX[0][i][0];
                        }

                        for (i = parseInt(ptsSX[1].length / 2); i < ptsSX[1].length; i++) {
                            sortedPtsArray[curvePart].push(ptsArray[ptsSX[1][i][0]][ptsSX[1][i][1]]);
                            lastPoint = ptsArray[ptsSX[1][i][0]][ptsSX[1][i][1]];
                            indexOfLastPoint = ptsSX[1][i][0];
                        }

                        for (i = parseInt(ptsSX[1].length / 2); i < ptsSX[1].length; i++) {
                            ptsArray[ptsSX[1][i][0]].splice(ptsSX[1][i][1], 1);
                        }

                        for (i = 0; i < ptsSX[0].length / 2; i++) {
                            ptsArray[ptsSX[0][i][0]].splice(ptsSX[0][i][1], 1);
                        }

                       	next = [];
                        if (ptsSX[1][ptsSX[1].length - 1] != undefined) {
                            for (i = 0; i < ptsArray[ptsSX[1][ptsSX[1].length - 1][0]].length; i++) {
                                if(equal(ptsArray[ptsSX[1][ptsSX[1].length-1][0]][i][1], lastPoint[1]+sY)) {
                                    next = [ptsSX[1][ptsSX[1].length - 1][0], i];
                                }
                            }
                        }

                        if (next.length != 0) {
                            sortedPtsArray[curvePart].push(ptsArray[next[0]][next[1]]);
                            lastPoint = ptsArray[next[0]][next[1]];
                            indexOfLastPoint = next[0];
                            ptsArray[next[0]].splice(next[1],1);
                        }
                    }
                }
            }

            // connect parts
            result = [];
            parts = 0;

            result[parts] = [];
            for (i = 0; i < sortedPtsArray[parts].length; i++) {
                result[parts].push(sortedPtsArray[parts][i]);
            }
            sortedPtsArray.splice(0, 1);

            while(sortedPtsArray.length > 0) {
                minDistance = 10000;
                idx = -1;
                isReverse = false;
				back = true;
                for (i = 0; i < sortedPtsArray.length; i++) {
                    if (result[parts][result[parts].length-1][0] > sortedPtsArray[i][0]) {
                         continue;
                    }

                    d = Geometry.distance(result[parts][result[parts].length-1],
                            sortedPtsArray[i][0], 2);
                    if (d < minDistance) {
                        minDistance = d;
                        idx = i;
                        back = true;
                        isReverse = false;
                    }

                    d = Geometry.distance(result[parts][result[parts].length-1],
                            sortedPtsArray[i][sortedPtsArray[i].length - 1], 2);
                    if (d < minDistance) {
                        minDistance = d;
                        idx = i;
                        back = true;
                        isReverse = true;
                    }

                    d = Geometry.distance(result[parts][0], sortedPtsArray[i][0], 2);
                    if (d < minDistance) {
                        minDistance = d;
                        idx = i;
                        back = false;
                        isReverse = false;
                    }

                    d = Geometry.distance(result[parts][0],
                        sortedPtsArray[i][sortedPtsArray[i].length - 1], 2);
                    if (d < minDistance) {
                        minDistance = d;
                        idx = i;
                        back = false;
                        isReverse = true;
                    }
                }

                if (minDistance > 4*stepX) {
                    result[++parts] = [];
                    idx = 0;
                }

                if (isReverse) {
                    sortedPtsArray[idx].reverse();
                }

                if (back) {
                    for (i = 0; i < sortedPtsArray[idx].length; i++) {
                        result[parts].push(sortedPtsArray[idx][i]);
                    }
                } else {
                    for (i=0; i < sortedPtsArray[idx].length; i++) {
                        result[parts].splice(0, 0, sortedPtsArray[idx][i]);
                    }
                }

                sortedPtsArray.splice(idx, 1);
            }

            return result;
        }
    };

    return Mat.Implicitplot;
};
