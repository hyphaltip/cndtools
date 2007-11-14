# Copyright (c) 2006
# Colin Dewey (University of Wisconsin-Madison)
# cdewey@biostat.wisc.edu
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

specialchars = {
    '"': "&quot;",
    '<': "&lt;",
    '>': "&gt;",
    '  ': " &nbsp;",
    '\n': "<br>",
    }

def escape(string):
    if isinstance(string, HTMLText):
        return str(string)
    else:
        string = str(string).replace("&", "&amp;")
        for special, replacement in specialchars.items():
            string = string.replace(special, replacement)
        return string

def _quote(string):
    return '"%s"' % string

def _makeattribute((name, value)):
    return "%s=%s" % (name, _quote(escape(value)))

def _makeattributes(attributes):
    if attributes:
        return ' ' + ' '.join(map(_makeattribute, attributes.items()))
    else:
        return ''

class HTMLText:
    def __init__(self, tag, singleton, noescape, *inside, **attributes):
        self.tag = tag
        self.singleton = singleton
        self.noescape = noescape
        self.inside = inside
        self.attributes = attributes
        
    def __str__(self):
        if self.singleton:
            return "<%s%s>" % (self.tag, _makeattributes(self.attributes))
        elif self.noescape:
            return ''.join(map(str, self.inside))
        else:
            text = ''.join(map(escape, self.inside))
            
            return "<%s%s>%s</%s>" % (self.tag,
                                      _makeattributes(self.attributes),
                                      text,
                                      self.tag)

elements = [
    ("a", 0),
    ("abbr", 0),
    ("acronym", 0),
    ("address", 0),
    ("applet", 0),
    ("area", 1),
    ("b", 0),
    ("base", 1),
    ("basefont", 1),
    ("bdo", 0),
    ("big", 0),
    ("blockquote", 0),
    ("body", 0),
    ("br", 1),
    ("button", 0),
    ("caption", 0),
    ("center", 0),
    ("cite", 0),
    ("code", 0),
    ("col", 1),
    ("colgroup", 0),
    ("dd", 0),
    ("html_del", 0),
    ("dfn", 0),
    ("dir", 0),
    ("div", 0),
    ("dl", 0),
    ("dt", 0),
    ("em", 0),
    ("fieldset", 0),
    ("font", 0),
    ("form", 0),
    ("frame", 1),
    ("frameset", 0),
    ("h1", 0),
    ("h2", 0),
    ("h3", 0),
    ("h4", 0),
    ("h5", 0),
    ("h6", 0),
    ("head", 0),
    ("hr", 1),
    ("html", 0),
    ("i", 0),
    ("iframe", 0),
    ("img", 1),
    ("input", 1),
    ("ins", 0),
    ("isindex", 1),
    ("kbd", 0),
    ("label", 0),
    ("legend", 0),
    ("li", 0),
    ("link", 1),
    ("html_map", 0),
    ("menu", 0),
    ("meta", 1),
    ("noframes", 0),
    ("noscript", 0),
    ("object", 0),
    ("ol", 0),
    ("optgroup", 0),
    ("option", 0),
    ("p", 0),
    ("param", 1),
    ("pre", 0),
    ("q", 0),
    ("s", 0),
    ("samp", 0),
    ("script", 0),
    ("select", 0),
    ("small", 0),
    ("span", 0),
    ("strike", 0),
    ("strong", 0),
    ("style", 0),
    ("sub", 0),
    ("sup", 0),
    ("table", 0),
    ("tbody", 0),
    ("td", 0),
    ("textarea", 0),
    ("tfoot", 0),
    ("th", 0),
    ("thead", 0),
    ("title", 0),
    ("tr", 0),
    ("tt", 0),
    ("u", 0),
    ("ul", 0),
    ("var", 0),
    ]

funcdef = """def %(tag)s(*inside, **attributes): return HTMLText('%(tag)s', %(singleton)d, 0, *inside, **attributes)"""

for tag, singleton in elements:
    exec(funcdef % vars())

def noescape(*inside):
    return HTMLText('', 0, 1, *inside)
