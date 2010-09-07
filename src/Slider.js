/*
    Copyright 2008,2009
        Matthias Ehmann,
        Michael Gerhaeuser,
        Carsten Miller,
        Bianca Valentin,
        Alfred Wassermann,
        Peter Wilfahrt

    This file is part of JSXGraph.

    JSXGraph is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    JSXGraph is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with JSXGraph.  If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * @fileoverview The geometry object Line is defined in this file. Line stores all
 * style and functional properties that are required to draw and move a line on
 * a board.
 */

/**
 * The Slider class is....
 * Slider (Schieberegler)
 * input: 3 arrays:
 * [x0,y0],[x1,y1],[min,start,max]
 * The slider is line from [x0,y0] to [x1,y1].
 * The position [x0,y0]  corresponds to the value "min",
 * [x1,y1] corresponds to the value max.
 * Initally, the slider is at position [x0,y0] + ([x1,y1]-[x0,y0])*start/(max-min)
 * The return value is an invisible point, whos X() or Y() value
 * returns the position between max and min,
 * Further, there is a method Value() returning the same value.
 * @class Creates a new basic slider object. Do not use this constructor to create a slider. Use {@link JXG.Board#create} with
 * type {@link Line}, {@link Arrow}, or {@link Axis} instead.  
 * Attributes: withTicks;
 * @constructor
 * @augments JXG.GeometryElement
 * @param {String,JXG.Board} board The board the new line is drawn on.
 * @param {Point} p1 Startpoint of the line.
 * @param {Point} p2 Endpoint of the line.
 * @param {String} id Unique identifier for this object. If null or an empty string is given,
 * an unique id will be generated by Board
 * @param {String} name Not necessarily unique name. If null or an
 * empty string is given, an unique name will be generated.
 * @see JXG.Board#generateName
 */
JXG.createSlider = function(board, parents, atts) {
    var pos0, pos1, smin, start, smax, sdiff, p1, p2, l1, ticks, ti, startx, starty, p3, l2, n, t,
        snapWidth, fixed;
        
    pos0 = parents[0];
    pos1 = parents[1];
    smin = parents[2][0];
    start = parents[2][1];
    smax = parents[2][2];
    sdiff = smax -smin;
    
    atts = JXG.checkAttributes(atts,{strokeColor:'#000000', fillColor:'#ffffff', withTicks:true});

    fixed = JXG.str2Bool(atts['fixed']);
    p1 = board.create('point', pos0, 
        {visible:!fixed, fixed:fixed, name:'',withLabel:false,face:'<>', size:5, strokeColor:'#000000', fillColor:'#ffffff'}); 
    p2 = board.create('point', pos1, 
        {visible:!fixed, fixed:fixed, name:'',withLabel:false,face:'<>', size:5, strokeColor:'#000000', fillColor:'#ffffff'}); 
    board.create('group',[p1,p2]);
    l1 = board.create('segment', [p1,p2], 
                {strokewidth:1, 
                name:'',
                withLabel:false,
                strokeColor:atts['strokeColor']});
    if (atts['withTicks']) {
        ticks  = 2;
        ti = board.create('ticks', [l1, p2.Dist(p1)/ticks],
                    {insertTicks:true, minorTicks:0, drawLabels:false, drawZero:true}); 
    }

    if (fixed) {
        p1.setProperty({needsRegularUpdate : false});
        p2.setProperty({needsRegularUpdate : false});
        l1.setProperty({needsRegularUpdate : false});
    }
    
    startx = pos0[0]+(pos1[0]-pos0[0])*(start-smin)/(smax-smin);
    starty = pos0[1]+(pos1[1]-pos0[1])*(start-smin)/(smax-smin);

    if (atts['snapWidth']!=null) snapWidth = atts['snapWidth'];
    if (atts['snapwidth']!=null) snapWidth = atts['snapwidth'];
    
    p3 = board.create('glider', [startx,starty,l1],
                {style:6,strokeColor:atts['strokeColor'],
                 fillColor:atts['fillColor'],
                 showInfobox:false,name:atts['name'], withLabel:false,
                 snapWidth:snapWidth});
    
    l2 = board.create('line', [p1,p3], 
                {straightFirst:false, 
                 straightLast:false, strokewidth:3, 
                 strokeColor:atts['strokeColor'],
                 name:'',
                 withLabel:false}); 
                 
    //p3.Value = function() { return this.position*(smax - smin)+smin; };
    //p3.type = JXG.OBJECT_TYPE_SLIDER;
    p3.Value = function() { return this.position*sdiff+smin; };
    p3._smax = smax;
    p3._smin = smin;

    if (typeof atts['withLabel']=='undefined' || atts['withLabel']==true) {
        if (atts['name'] && atts['name']!='') {
            n = atts['name'] + ' = ';
        } else {
            n = '';
        }
        t = board.create('text', [function(){return (p2.X()-p1.X())*0.05+p2.X();},
                                  function(){return (p2.Y()-p1.Y())*0.05+p2.Y();},
                                  function(){return n+(p3.Value()).toFixed(2);}],
                                     {name:''}); 
    }                                     
    return p3;
};    

JXG.JSXGraph.registerElement('slider', JXG.createSlider);
