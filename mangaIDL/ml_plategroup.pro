;+
; NAME:
;   ml_plategroup
;
; PURPOSE:
;   Returns the plate group name for a given plate.
;   This groups a set of 100 plate IDs together for use in keeping
;   a directory structure clean.
;   E.g., the group name for plate 6534 would be '0065XX'.
;   The group name for plate 12432 would be '0124XX'.
;
; CALLING SEQUENCE:
;   groupname=ml_plategroup(plateid)
;
; INPUTS:
;   plateid - long integer plate ID
;
; OPTIONAL INPUTS:
;   /design - Input integer is design ID, not plate ID.
;             In this case, output is (e.g.) 'D65XX'
;
; OUTPUTS:
;   groupname- String format plate group name.
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   October-2013  Written by David Law (Dunlap Institute; drlaw@di.utoronto.ca)
;   03-Feb-2014 Add leading '00' for plateid<10k for consistency with SDSS III
;-

function ml_plategroup,plateid,design=design

; Convert to string format plateid
groupname=strcompress(string(plateid),/remove_all)

; We don't have any plate IDs less than 1000, so don't
; worry about them.  Just replace end 2 characters with 'X'
strput,groupname,'X',strlen(groupname)-1
strput,groupname,'X',strlen(groupname)-2

; Buffer with leading zeroes
if (long(plateid) le 9999) then groupname='00'+groupname
if (long(plateid) ge 10000) then groupname='0'+groupname

if (keyword_set(design)) then groupname='D'+groupname

return,groupname
end
