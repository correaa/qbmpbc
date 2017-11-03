\begingroup
\obeyspaces\obeylines\global\let^^M=\jsR%
\catcode`\"=12%
\gdef\dljspgfplotsJSiii{%
/*********************************************************************************
 * function sprintf() - written by Kevin van Zonneveld as part of the php to javascript
 * conversion project.
 *
 * More info at: http://kevin.vanzonneveld.net/techblog/article/phpjs_licensing/
 *
 * This is version: 1.33
 * php.js is copyright 2008 Kevin van Zonneveld.
 *
 * Portions copyright Michael White (http://crestidg.com), _argos, Jonas
 * Raoni Soares Silva (http://www.jsfromhell.com), Legaev Andrey, Ates Goral
 * (http://magnetiq.com), Philip Peterson, Martijn Wieringa, Webtoolkit.info
 * (http://www.webtoolkit.info/), Carlos R. L. Rodrigues
 * (http://www.jsfromhell.com), Ash Searle (http://hexmen.com/blog/),
 * Erkekjetter, GeekFG (http://geekfg.blogspot.com), Johnny Mast
 * (http://www.phpvrouwen.nl), marrtins, Alfonso Jimenez
 * (http://www.alfonsojimenez.com), Aman Gupta, Arpad Ray
 * (mailto:arpad@php.net), Karol Kowalski, Mirek Slugen, Thunder.m, Tyler
 * Akins (http://rumkin.com), d3x, mdsjack (http://www.mdsjack.bo.it), Alex,
 * Alexander Ermolaev (http://snippets.dzone.com/user/AlexanderErmolaev),
 * Allan Jensen (http://www.winternet.no), Andrea Giammarchi
 * (http://webreflection.blogspot.com), Arno, Bayron Guevara, Ben Bryan,
 * Benjamin Lupton, Brad Touesnard, Brett Zamir, Cagri Ekin, Cord, David,
 * David James, DxGx, FGFEmperor, Felix Geisendoerfer
 * (http://www.debuggable.com/felix), FremyCompany, Gabriel Paderni, Howard
 * Yeend, J A R, Leslie Hoare, Lincoln Ramsay, Luke Godfrey, MeEtc
 * (http://yass.meetcweb.com), Mick@el, Nathan, Nick Callen, Ozh, Pedro Tainha
 * (http://www.pedrotainha.com), Peter-Paul Koch
 * (http://www.quirksmode.org/js/beat.html), Philippe Baumann, Sakimori,
 * Sanjoy Roy, Simon Willison (http://simonwillison.net), Steve Clay, Steve
 * Hilder, Steven Levithan (http://blog.stevenlevithan.com), T0bsn, Thiago
 * Mata (http://thiagomata.blog.com), Tim Wiel, XoraX (http://www.xorax.info),
 * Yannoo, baris ozdil, booeyOH, djmix, dptr1988, duncan, echo is bad, gabriel
 * paderni, ger, gorthaur, jakes, john (http://www.jd-tech.net), kenneth,
 * loonquawl, penutbutterjelly, stensi
 *
 * Dual licensed under the MIT (MIT-LICENSE.txt)
 * and GPL (GPL-LICENSE.txt) licenses.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL KEVIN VAN ZONNEVELD BE LIABLE FOR ANY CLAIM, DAMAGES
 * OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
// ATTENTION: this method has been masked such that special characters of TeX and javascript
// don't produce problems.
function sprintf( ) {
    // Return a formatted string
    //
    // +    discuss at: http://kevin.vanzonneveld.net/techblog/article/javascript_equivalent_for_phps_sprintf/
    // +       version: 804.1712
    // +   original by: Ash Searle (http://hexmen.com/blog/)
    // + namespaced by: Michael White (http://crestidg.com)
    // *     example 1: sprintf("\pgfplotsPERCENT01.2f", 123.1);
    // *     returns 1: 123.10

    var regex = /\pgfplotsPERCENT\pgfplotsPERCENT\pgfplotsVERTBAR\pgfplotsPERCENT(\d+\$)?([-+\pgfplotsHASH0 ]*)(\*\d+\$\pgfplotsVERTBAR\*\pgfplotsVERTBAR\d+)?(\.(\*\d+\$\pgfplotsVERTBAR\*\pgfplotsVERTBAR\d+))?([scboxXuidfegEG])/g;
    var a = arguments, i = 0, format = a[i++];

    // pad()
    var pad = function(str, len, chr, leftJustify) {
        var padding = (str.length >= len) ? '' : Array(1 + len - str.length >>> 0).join(chr);
        return leftJustify ? str + padding : padding + str;
    };

    // justify()
    var justify = function(value, prefix, leftJustify, minWidth, zeroPad) {
        var diff = minWidth - value.length;
        if (diff > 0) {
            if (leftJustify \pgfplotsVERTBAR\pgfplotsVERTBAR !zeroPad) {
            value = pad(value, minWidth, ' ', leftJustify);
            } else {
            value = value.slice(0, prefix.length) + pad('', diff, '0', true) + value.slice(prefix.length);
            }
        }
        return value;
    };

    // formatBaseX()
    var formatBaseX = function(value, base, prefix, leftJustify, minWidth, precision, zeroPad) {
        // Note: casts negative numbers to positive ones
        var number = value >>> 0;
        prefix = prefix && number && {'2': '0b', '8': '0', '16': '0x'}[base] \pgfplotsVERTBAR\pgfplotsVERTBAR '';
        value = prefix + pad(number.toString(base), precision \pgfplotsVERTBAR\pgfplotsVERTBAR 0, '0', false);
        return justify(value, prefix, leftJustify, minWidth, zeroPad);
    };

    // formatString()
    var formatString = function(value, leftJustify, minWidth, precision, zeroPad) {
        if (precision != null) {
            value = value.slice(0, precision);
        }
        return justify(value, '', leftJustify, minWidth, zeroPad);
    };

    // finalFormat()
    var doFormat = function(substring, valueIndex, flags, minWidth, _, precision, type) {
        if (substring == '\pgfplotsPERCENT\pgfplotsPERCENT') return '\pgfplotsPERCENT';

        // parse flags
        var leftJustify = false, positivePrefix = '', zeroPad = false, prefixBaseX = false;
        for (var j = 0; flags && j < flags.length; j++) switch (flags.charAt(j)) {
            case ' ': positivePrefix = ' '; break;
            case '+': positivePrefix = '+'; break;
            case '-': leftJustify = true; break;
            case '0': zeroPad = true; break;
            case '\pgfplotsHASH': prefixBaseX = true; break;
        }

        // parameters may be null, undefined, empty-string or real valued
        // we want to ignore null, undefined and empty-string values
        if (!minWidth) {
            minWidth = 0;
        } else if (minWidth == '*') {
            minWidth = +a[i++];
        } else if (minWidth.charAt(0) == '*') {
            minWidth = +a[minWidth.slice(1, -1)];
        } else {
            minWidth = +minWidth;
        }

        // Note: undocumented perl feature:
        if (minWidth < 0) {
            minWidth = -minWidth;
            leftJustify = true;
        }

        if (!isFinite(minWidth)) {
            throw new Error('sprintf: (minimum-)width must be finite');
        }

        if (!precision) {
            precision = 'fFeE'.indexOf(type) > -1 ? 6 : (type == 'd') ? 0 : void(0);
        } else if (precision == '*') {
            precision = +a[i++];
        } else if (precision.charAt(0) == '*') {
            precision = +a[precision.slice(1, -1)];
        } else {
            precision = +precision;
        }

        // grab value using valueIndex if required?
        var value = valueIndex ? a[valueIndex.slice(0, -1)] : a[i++];

        switch (type) {
            case 's': return formatString(String(value), leftJustify, minWidth, precision, zeroPad);
            case 'c': return formatString(String.fromCharCode(+value), leftJustify, minWidth, precision, zeroPad);
            case 'b': return formatBaseX(value, 2, prefixBaseX, leftJustify, minWidth, precision, zeroPad);
            case 'o': return formatBaseX(value, 8, prefixBaseX, leftJustify, minWidth, precision, zeroPad);
            case 'x': return formatBaseX(value, 16, prefixBaseX, leftJustify, minWidth, precision, zeroPad);
            case 'X': return formatBaseX(value, 16, prefixBaseX, leftJustify, minWidth, precision, zeroPad).toUpperCase();
            case 'u': return formatBaseX(value, 10, prefixBaseX, leftJustify, minWidth, precision, zeroPad);
            case 'i':
            case 'd': {
                        var number = parseInt(+value);
                        var prefix = number < 0 ? '-' : positivePrefix;
                        value = prefix + pad(String(Math.abs(number)), precision, '0', false);
                        return justify(value, prefix, leftJustify, minWidth, zeroPad);
                    }
            case 'e':
            case 'E':
            case 'f':
            case 'F':
            case 'g':
            case 'G':
                        {
                        var number = +value;
                        var prefix = number < 0 ? '-' : positivePrefix;
                        var method = ['toExponential', 'toFixed', 'toPrecision']['efg'.indexOf(type.toLowerCase())];
                        var textTransform = ['toString', 'toUpperCase']['eEfFgG'.indexOf(type) \pgfplotsPERCENT 2];
                        value = prefix + Math.abs(number)[method](precision);
                        return justify(value, prefix, leftJustify, minWidth, zeroPad)[textTransform]();
                    }
            default: return substring;
        }
    };

    return format.replace(regex, doFormat);
}
/*********************************************************************************/


var lastMouseX = -1;
var lastMouseY = -1;
var posOnMouseDownX = -1;
var posOnMouseDownY = -1;

// preallocation.
var tmpArray1 = new Array(2);
var tmpArray2 = new Array(2);

/**
 * Takes an already existing TextField, changes its value to \c value and places it at (x,y).
 * Additional \c displayOpts will be used to format it.
 */
function initTextField( x,y, textField, displayOpts, value )
{
textField.value = ""+value;
var R = textField.rect;
R[0] = x;
R[1] = y;
R[2] = R[0] + displayOpts.textSize/2*Math.max( 5,value.length );
R[3] = R[1] - 1.5*displayOpts.textSize;
textField.rect = R;
textField.textFont = displayOpts.textFont;
textField.textSize = displayOpts.textSize;
textField.fillColor = displayOpts.fillColor;//["RGB",1,1,.855];
textField.doNotSpellCheck = true;
textField.readonly = true;
if( displayOpts.printable )
textField.display = display.visible;
else
textField.display = display.noPrint;
}

/**
 * Changes all required Field values of \c plotRegionField, inserts the proper
 * value and displays it at the pdf positions (x,y) .
 *
 * @param plotRegionField a reference to a Field object.
 * @param x the x coordinate where the annotation shall be placed and which is used to determine
 *  the annotation text.
 * @param y the corresponding y coord.
 * @param axisAnnotObj An object containing axis references.
 * @param displayOpts An object for display flags.
 * @param[out] retCoords will be filled with the point in axis coordinates.
 */
function placeAnnot( plotRegionField, x,y, textField, axisAnnotObj, displayOpts, retCoords )
{
// Get and modify bounding box. The mouse movement is only accurate up to one point
// (mouseX and mouseY are integers), so the bounding box should be an integer as well.
var R = plotRegionField.rect;
R[0] = Math.round(R[0]);
R[1] = Math.round(R[1]);
R[2] = Math.round(R[2]);
R[3] = Math.round(R[3]);
plotRegionField.rect= R;

var w = R[2] - R[0];
var h = R[1] - R[3];

// compute position 0 <= mx <= 1, 0<= my <= 1 relative to lower(!) left corner.
var mx = x - R[0];
var my = h - (R[1] - y);

var unitx = mx / w;
var unity = my / h;

var realx = axisAnnotObj.xmin + unitx * (axisAnnotObj.xmax - axisAnnotObj.xmin);
var realy = axisAnnotObj.ymin + unity * (axisAnnotObj.ymax - axisAnnotObj.ymin);

var transformedCoordx = realx;
var transformedCoordy = realy;

if( axisAnnotObj.xscale.length >= 3 && axisAnnotObj.xscale.substr(0,3) == "log" ) {
if( axisAnnotObj.xscale.length > 4 ) // log:<basis>
realx = realx * Math.log( axisAnnotObj.xscale.substr(4) );
else {
// pgfplots handles log plots base e INTERNALLY, but uses base 10 for display.
// convert to base 10:
transformedCoordx = realx / Math.log(10);
}
realx = Math.exp(realx);
}
if( axisAnnotObj.yscale.length >= 3 && axisAnnotObj.yscale.substr(0,3) == "log" ) {
if( axisAnnotObj.yscale.length > 4 ) // log:<basis>
realy = realy * Math.log( axisAnnotObj.yscale.substr(4) );
else {
// pgfplots handles log plots base e INTERNALLY, but uses base 10 for display.
// convert to base 10:
transformedCoordy = realy / Math.log(10);
}
realy = Math.exp(realy);
}

// console.println( "w = " + w + "; h = " + h );
// console.println( "mx = " + mx + "; my = " + my );
// console.println( "unitx = " + unitx + "; unity " + unity );

initTextField( x,y, textField, displayOpts,
//util.printf( "(\pgfplotsPERCENT.2f,\pgfplotsPERCENT.2f)", realx,realy )
sprintf( displayOpts.pointFormat, realx,realy) );

if( retCoords ) {
retCoords[0] = transformedCoordx;
retCoords[1] = transformedCoordy;
}

}

/**
 * @param formName the name of the clickable button. It is expected to be as large as the underlying plot.
 * @param axisAnnotObj an object with the fields
 *   - xmin, xmax
 *   - ymin, ymax
 *   - xscale, yscale
 * @param displayOpts an object with the fields
 *   - pointFormat an sprintf format string to format the final point coordinates.
 *   The default is  "(\pgfplotsPERCENT.2f,\pgfplotsPERCENT.2f)"
 *   - fillColor the fill color for the annotation. Options are
 *    transparent, gray, RGB or CMYK color. Default is
 *       ["RGB",1,1,.855]
 *  - textFont / textSize
 */
function processAnnotatedPlot(formName, axisAnnotObj, displayOpts, bMouseUp )
{
if( !bMouseUp ) {
posOnMouseDownX = mouseX;
posOnMouseDownY = mouseY;
return;
}
var result = this.getField( formName + "-result");
var result2 = this.getField( formName + "-result2");
var slope  = this.getField( formName + "-slope" );
if( !result ) {
console.println( "WARNING: there is no TextField \"" + formName + "-result\" to display results for interactive element \"" + formName + "\"");
return;
}
result2.display = display.hidden;
slope.display = display.hidden;

// clicking twice onto the same point hides it:
if( result.display != display.hidden &&
Math.abs(lastMouseX - mouseX) < 5 &&
Math.abs(lastMouseY - mouseY) < 5 )
{
result.display = display.hidden;
return;
}
lastMouseX = mouseX;
lastMouseY = mouseY;

var a = this.getField(formName);
if( ! a ) {
console.println( "Warning: there is no form named \"" + formName + "\"" );
return;
}
if( Math.abs( lastMouseX - posOnMouseDownX ) > 6 \pgfplotsVERTBAR\pgfplotsVERTBAR
Math.abs( lastMouseY - posOnMouseDownY ) > 6 )
{
// dragging the mouse results in slope computation:
// placeAnnot shows the endpoint coords and returns the (transformed) coordinates into tmpArray1 and tmpArray2:
placeAnnot( a, posOnMouseDownX, posOnMouseDownY, result, axisAnnotObj, displayOpts, tmpArray1 );
placeAnnot( a, lastMouseX, lastMouseY, result2, axisAnnotObj, displayOpts, tmpArray2 );

var m =  ( tmpArray2[1] - tmpArray1[1] ) / ( tmpArray2[0] - tmpArray1[0] );
var n =  tmpArray1[1] - m * tmpArray1[0];

initTextField(
0.5 * ( lastMouseX + posOnMouseDownX ),
0.5 * ( lastMouseY + posOnMouseDownY ),
slope,
displayOpts,
sprintf( displayOpts.slopeFormat, m, n ) );

// FIXME! these document rights seem to forbid modifications to annotations, although they work for text fields!?
//var lineobj = this.getAnnot( a.page, formName + '-line' );
//console.println( 'lineobj = ' + lineobj );
//lineobj.points = [[lastMouseX,lastMouseY],[posOnMouseDownX,posOnMouseDownY]];
//lineobj.display = display.visible;

} else {
placeAnnot( a, lastMouseX, lastMouseY, result, axisAnnotObj, displayOpts, null );
}
}
}%
\endgroup
\begingroup 
\ccpdftex%
\input{dljscc.def}%
\immediate\pdfobj{ << /S /JavaScript /JS (\dljspgfplotsJSiii) >> }
\xdef\objpgfplotsJSiii{\the\pdflastobj\space0 R}
\endgroup 
\endinput
