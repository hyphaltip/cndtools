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

import telnetlib
import re

class TimeoutError(Exception):
    pass

def _waitForPrompt(tn, prompt, timeout):
    """Wait for string PROMPT on telnet connection TN until
       TIMEOUT is reached

       Raises TimeoutError if prompt is not read before TIMEOUT
    """
    
    s = tn.read_until(prompt, timeout)
    if s.find(prompt) == -1:
        raise TimeoutError, "GLADIS did not respond in %d sec" % timeout
    return s

def renewAllBooks(SID, SSN, timeout=10):
    """Renew all books checked out by user with specified SID and SSN"""
    
    HOST = "gladis.berkeley.edu" # GLADIS server
    PROMPT = "===> " # Standard GLADIS prompt
    RETURN = "\r"
    ITEMS_PER_SCREEN = 7 # Number of books shown on inventory screen at a time

    tn = telnetlib.Telnet(HOST)

    # Log into GLADIS
    s = _waitForPrompt(tn, "Enter Choice> ", timeout)
    tn.write("GLADIS" + RETURN)
    s = _waitForPrompt(tn, PROMPT, timeout)
    tn.write(RETURN)
    s = _waitForPrompt(tn, PROMPT, timeout)
    
    # Enter inventory section with SID and SSN
    tn.write("inv" + RETURN)
    s = _waitForPrompt(tn, PROMPT, timeout)
    tn.write(str(SID) + RETURN)
    s = _waitForPrompt(tn, PROMPT, timeout)
    tn.write(str(SSN) + RETURN)
    s = _waitForPrompt(tn, PROMPT, timeout)
    
    # Read number of books checkout from inventory screen        
    m = re.search(r"(\d+) items? checked out\.", s)
    items = int(m.group(1))

    renewPat = re.compile(r"\*\*>  (Item # \d+ SUCCESSFULLY RENEWED\.  " +
                          r"It is now due \d+/\d+/\d+\.)")
    results = ["Unsuccessful"] * items

    # Renew all books    
    for i in xrange(1, items + 1):
        # Jump to next screen if done with books on current screen
        if (i != 1) and (i % ITEMS_PER_SCREEN == 1):
            tn.write(RETURN)
            s = _waitForPrompt(tn, PROMPT, timeout)

        # Renew book i            
        tn.write(("ren %d" + RETURN) % i)
        s = _waitForPrompt(tn, PROMPT, timeout)
        m = renewPat.search(s)
        if m is not None:
            results[i - 1] = m.group(1)

    # Log off of GLADIS
    tn.write("ST" + RETURN)
    tn.close()

    return results    
