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
  //base/constants
  //math/math
  //math/geometry
  //math/numerics
  //utils/type
   elements:
    //point
    curve
  */

 /**
  * @fileoverview In this file the conic sections defined.
  */

define([
     'jxg', 'base/constants', 'base/coords', 'math/math', 'math/numerics', 'math/geometry', 'utils/type', 'base/point', 'base/curve'
], function (JXG, Const, Coords, Mat, Numerics, Geometry, Type, Point, Curve) {

    JXG.createImplicitPlot = function (board, parents, attributes) {
        var el, attr;

        attr = JXG.copyAttributes(attributes, board.options, 'curve');
        el = board.create('curve', [[0], [0]], attr);

        el.equation = parents[0];
        el.updateDataArray = function() {
            var x, y, i, j, len,
                result,
                val = this.equation(0, 0),
                bb = this.board.getBoundingBox(),
                array = [],
                searchDepth = 8,
                d = searchDepth + 2,
                step = (bb[2] - bb[0]) / (d * d);

            this.dataX = [];
            this.dataY = [];

            find_points(this.equation, array, 0, bb[0], bb[3], bb[2]-bb[0], searchDepth);
            result = organizePoints(this.equation, array, step);
            len = result.length;
            console.log(len);

            for (i = 0; i < len; i++) {
                for (j = 0; j < result[i].length; j++) {
                    this.dataX.push(result[i][j][0]);
                    this.dataY.push(result[i][j][1]);
                }
                this.dataX.push(NaN);
                this.dataY.push(NaN);
            }
        };

        el.prepareUpdate().update().updateRenderer();
        return el;
    };

    JXG.registerElement('implicitplot', JXG.createImplicitPlot);

    return {
        createImplicitPlot: JXG.createImplicitPlot
    };

});


function equal(v1, v2)
{
    var EPS = 10e-6;
    return Math.abs(v1-v2) < EPS;
}

function addInArray(ptsArray, point)
{
    if (ptsArray.length == 0)
    {
        ptsArray[0] = new Array();
        ptsArray[0].push(point);
    }
    else
    {
        var pos;
        for(pos=0; pos < ptsArray.length && ptsArray[pos][0][0] < point[0]; pos++)
            ;
        if (pos == ptsArray.length)
            ptsArray[pos] = new Array();
        else if(point[0] != ptsArray[pos][0][0])
            ptsArray.splice(pos, 0, new Array());
        ptsArray[pos].push(point);
    }
}

function distance(p1, p2)
{
    return Math.sqrt(Math.pow(p1[0]-p2[0], 2) + Math.pow(p1[1]-p2[1], 2));
}

function organizePoints(f, ptsArray, step)
{
    if (ptsArray.length == 0)
        return [];

    sortedPtsArray = new Array();
    curvePart = 0;
    // find curve parts

    var stepX = step;
    var stepY = step;

    sortedPtsArray[curvePart] = new Array();
    sortedPtsArray[curvePart].push(ptsArray[0][0]);
    var lastPoint = ptsArray[0][0];
    var indexOfLastPoint = 0;

    ptsArray[0].splice(0,1);
    if (ptsArray[0].length == 0)
        ptsArray.splice(0,1);

    var noPoints = 0;
    for(var i=0; i<ptsArray.length; i++)
        for(var j=0; j<ptsArray[i].length; j++)
            noPoints++;

    for(var iter=0; iter<noPoints; iter++)
    {
        var neighbours = new Array();

        var previousIndex =
                (indexOfLastPoint != 0 && ptsArray[indexOfLastPoint-1].length != 0 && equal(ptsArray[indexOfLastPoint-1][0][0] + stepX, lastPoint[0]) ) ?
                    indexOfLastPoint-1 :
                    -1;
        var nextIndex = (indexOfLastPoint != ptsArray.length-1 && ptsArray[indexOfLastPoint+1].length != 0 && equal(ptsArray[indexOfLastPoint+1][0][0] - stepX, lastPoint[0]) ) ?
                    indexOfLastPoint+1 :
                    -1;

        for(var i=0; i<ptsArray[indexOfLastPoint].length; i++)
            if (equal(lastPoint[1], ptsArray[indexOfLastPoint][i][1]+stepY) ||
                equal(lastPoint[1], ptsArray[indexOfLastPoint][i][1]-stepY) )
                neighbours.push([indexOfLastPoint, i]);

        if (previousIndex != -1)
            for(var i=0; i<ptsArray[previousIndex].length; i++)
                if (equal(lastPoint[1], ptsArray[previousIndex][i][1]))
                    neighbours.push([previousIndex, i]);

        if (nextIndex != -1)
            for(var i=0; i<ptsArray[nextIndex].length; i++)
                if (equal(lastPoint[1], ptsArray[nextIndex][i][1]))
                    neighbours.push([nextIndex, i]);

        if (neighbours.length == 0)
        {
            var pos = 0;
            for(pos=0; pos<ptsArray.length; pos++)
                if (ptsArray[pos].length != 0)
                    break;
            if (pos == ptsArray.length)
                break;

            var noPossiblePoints = ptsArray[pos].length,
                minYin = 0,
                maxYin = 0,
                startIndex = 0;

            if (noPossiblePoints > 1)
            {
                var neighbourIndex = 0;
                for(var i=1; i<noPossiblePoints; i++)
                    if (ptsArray[pos][i][1] < ptsArray[pos][minYin][1] )
                    {
                        minYin = i;
                        maxYin = i;
                    }
                for(var i=0; i<noPossiblePoints; i++)
                {
                    if (equal(ptsArray[pos][maxYin][1], ptsArray[pos][i][1]-stepY))
                        maxYin = i;
                }

                if (pos+1 < ptsArray.length-1 && ptsArray[pos+1].length > 1)
                {
                    for(var i=0; i<ptsArray[pos+1].length; i++)
                    {
                        if (equal(ptsArray[pos][minYin][1], ptsArray[pos+1][i][1]))
                        {
                            neighbourIndex = minYin;
                            break;
                        }
                        if (equal(ptsArray[pos][maxYin][1], ptsArray[pos+1][i][1]))
                        {
                            neighbourIndex = maxYin;
                            break;
                        }
                    }
                }
                startIndex = (neighbourIndex == maxYin) ? minYin : maxYin;
            }
            sortedPtsArray[++curvePart] = new Array();
            sortedPtsArray[curvePart].push(ptsArray[pos][startIndex]);

            lastPoint = ptsArray[pos][startIndex];
            indexOfLastPoint = pos;
            ptsArray[pos].splice(startIndex,1);
        }
        else
        {
            if (sortedPtsArray[curvePart].length == 1)
            { // start of closed curve
                sortedPtsArray[curvePart].push(ptsArray[neighbours[0][0]][neighbours[0][1]]);
                lastPoint = ptsArray[neighbours[0][0]][neighbours[0][1]];
                indexOfLastPoint = neighbours[0][0];
                ptsArray[neighbours[0][0]].splice(neighbours[0][1],1);
            }
            else if (neighbours.length == 1)
            { // just one neighbour point
                sortedPtsArray[curvePart].push(ptsArray[neighbours[0][0]][neighbours[0][1]]);
                lastPoint = ptsArray[neighbours[0][0]][neighbours[0][1]];
                indexOfLastPoint = neighbours[0][0];
                ptsArray[neighbours[0][0]].splice(neighbours[0][1],1);

                if (sortedPtsArray[curvePart].length > 3)
                {
                    var lastIndex = sortedPtsArray[curvePart].length-1;
                    if ( (  equal(sortedPtsArray[curvePart][lastIndex-2][0], sortedPtsArray[curvePart][lastIndex-1][0]) &&
                            equal(sortedPtsArray[curvePart][lastIndex-1][0], sortedPtsArray[curvePart][lastIndex][0] )) ||
                         (  equal(sortedPtsArray[curvePart][lastIndex-2][1], sortedPtsArray[curvePart][lastIndex-1][1]) &&
                            equal(sortedPtsArray[curvePart][lastIndex-1][1], sortedPtsArray[curvePart][lastIndex][1]) ) )
                                sortedPtsArray[curvePart].splice(lastIndex-1, 1);
                }
            }
            else if (neighbours.length > 1)
            {
                var noX = 0,
                    noY = 0,
                    sX, sY,
                    ptsSX = new Array();

                if (equal(lastPoint[1], ptsArray[neighbours[0][0]][neighbours[0][1]][1]))
                    sX = ptsArray[neighbours[0][0]][neighbours[0][1]][0] - lastPoint[0];
                else
                    sX = ptsArray[neighbours[1][0]][neighbours[1][1]][0] - lastPoint[0];
                if (equal(lastPoint[0], ptsArray[neighbours[0][0]][neighbours[0][1]][0]))
                    sY = ptsArray[neighbours[0][0]][neighbours[0][1]][1] - lastPoint[1];
                else
                    sY = ptsArray[neighbours[1][0]][neighbours[1][1]][1] - lastPoint[1];

                ptsSX[0] = new Array();
                ptsSX[1] = new Array();

                var signSX = (sX > 0) ? 1 : -1;


                for(var i=indexOfLastPoint+signSX; i<ptsArray.length && i>=0; i+=signSX)
                {
                    noX=0;
                    for(var j=0; j<ptsArray[i].length; j++)
                    {
                        if (equal(lastPoint[1], ptsArray[i][j][1]))
                        {
                            noX++;
                            ptsSX[0].push([i, j]);
                        }
                        if (equal(lastPoint[1]+sY, ptsArray[i][j][1]))
                        {
                            ptsSX[1].push([i, j]);
                        }
                    }
                    if (noX == 0)
                        break;
                }

                for(var i=0; i<ptsSX[0].length/2; i++)
                {
                    sortedPtsArray[curvePart].push(ptsArray[ptsSX[0][i][0]][ptsSX[0][i][1]]);
                    lastPoint = ptsArray[ptsSX[0][i][0]][ptsSX[0][i][1]];
                    indexOfLastPoint = ptsSX[0][i][0];
                }
                for(var i=parseInt(ptsSX[1].length/2); i<ptsSX[1].length; i++)
                {
                    sortedPtsArray[curvePart].push(ptsArray[ptsSX[1][i][0]][ptsSX[1][i][1]]);
                    lastPoint = ptsArray[ptsSX[1][i][0]][ptsSX[1][i][1]];
                    indexOfLastPoint = ptsSX[1][i][0];
                }
                for(var i=parseInt(ptsSX[1].length/2); i<ptsSX[1].length; i++)
                    ptsArray[ptsSX[1][i][0]].splice(ptsSX[1][i][1],1);
                for(var i=0; i<ptsSX[0].length/2; i++)
                    ptsArray[ptsSX[0][i][0]].splice(ptsSX[0][i][1],1);

                var next = new Array();
                if (ptsSX[1][ptsSX[1].length-1] != undefined)
                    for(var i=0; i<ptsArray[ptsSX[1][ptsSX[1].length-1][0]].length; i++)
                    {
                        if(equal(ptsArray[ptsSX[1][ptsSX[1].length-1][0]][i][1], lastPoint[1]+sY))
                            next = [ptsSX[1][ptsSX[1].length-1][0], i];
                    }

                if(next.length != 0)
                {
                    sortedPtsArray[curvePart].push(ptsArray[next[0]][next[1]]);
                    lastPoint = ptsArray[next[0]][next[1]];
                    indexOfLastPoint = next[0];
                    ptsArray[next[0]].splice(next[1],1);
                }
            }
        }
    }

    // connect parts
    var result = new Array();
    parts = 0;

    result[parts] = new Array();
    for(var i=0; i<sortedPtsArray[parts].length; i++)
        result[parts].push(sortedPtsArray[parts][i]);
    sortedPtsArray.splice(0, 1);

    while(sortedPtsArray.length > 0)
    {
        var minDistance = 10000;
        var idx = -1;
        var reverse = false;
        var back = true;
        for(i=0; i<sortedPtsArray.length; i++)
        {
            if(result[parts][result[parts].length-1][0] > sortedPtsArray[i][0])
                 continue;
            if (distance(result[parts][result[parts].length-1], sortedPtsArray[i][0]) < minDistance)
            {
                minDistance = distance(result[parts][result[parts].length-1], sortedPtsArray[i][0]);
                idx = i;
                back = true;
                reverse = false;
            }
            if (distance(result[parts][result[parts].length-1], sortedPtsArray[i][sortedPtsArray[i].length-1]) < minDistance)
            {
                minDistance = distance(result[parts][result[parts].length-1], sortedPtsArray[i][sortedPtsArray[i].length-1]);
                idx = i;
                back = true;
                reverse = true;
            }
            if (distance(result[parts][0], sortedPtsArray[i][0]) < minDistance)
            {
                minDistance = distance(result[parts][0], sortedPtsArray[i][0]);
                idx = i;
                back = false;
                reverse = false;
            }
            if (distance(result[parts][0], sortedPtsArray[i][sortedPtsArray[i].length-1]) < minDistance)
            {
                minDistance = distance(result[parts][0], sortedPtsArray[i][sortedPtsArray[i].length-1]);
                idx = i;
                back = false;
                reverse = true;
            }
        }
        if (minDistance > 4*stepX)
        {
            result[++parts] = new Array();
            idx = 0;
        }

        if (reverse)
            sortedPtsArray[idx].reverse();
        if (back)
            for (var i=0; i<sortedPtsArray[idx].length; i++)
                result[parts].push(sortedPtsArray[idx][i]);
        else
            for (var i=0; i<sortedPtsArray[idx].length; i++)
                result[parts].splice(0, 0, sortedPtsArray[idx][i]);

        sortedPtsArray.splice(idx, 1);
    }

    return result;
}

function find_points(f, array, depth, x, y, d, searchDepth)
{
    var firstDepth = searchDepth;
    var secondDepth = searchDepth+2;

    if (depth < firstDepth)
        subdivide(f, array, depth, x, y, d, searchDepth);
    else if (contour_present(f, x, y, d))
    {
        if (depth < secondDepth)
            subdivide(f, array, depth, x, y, d, searchDepth);
        else
        {
            var px = (2*x+d)/2;
            var py = (2*y+d)/2;
            point = [px, py];
            addInArray(array, point);
        }
    }
}

function subdivide(f, array, depth, x, y, d, searchDepth)
{
    find_points(f, array, depth+1, x, y, d/2, searchDepth);
    find_points(f, array, depth+1, x+d/2, y, d/2, searchDepth);
    find_points(f, array, depth+1, x+d/2, y+d/2, d/2, searchDepth);
    find_points(f, array, depth+1, x, y+d/2, d/2, searchDepth);
}

function contour_present(f, x, y, d)
{
    if (f(x,y) * f(x+d, y) < 0 || f(x, y) * f(x, y+d) < 0 || f(x+d, y) * f(x+d, y+d) < 0)
        return true;
    else if (f(x,y) == 0 && ( f(x+d,y) * f(x+d,y+d) < 0 || f(x,y+d) * f(x+d,y+d) < 0 || f(x+d, y+d) == 0) )
        return true;
    else if (f(x+d,y) == 0 && ( f(x,y) * f(x,y+d) < 0 || f(x,y+d) * f(x+d,y+d) < 0 || f(x, y+d) == 0) )
        return true;
    return false;
}
