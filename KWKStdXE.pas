unit KWKStdXE;
{(c) Copyright 1986-2012, Keith W. Kintigh}

interface

uses ansistrings, strutils;

  type   matrix_elmt=double;

const  missing=1.0e32;
       max_antana_errors=20;
       tokenlen=32;
       years='1986-2007';
       label_width=12;

       {database definition}
         maxvlabellen: integer=12;
         maxclabellen=12;
       {end of database definition constants}

type Int_Vector = array of integer;
     Real_Vector = array of real;
     String_Vector = array of string;
     Int_Array_2D = array of array of integer;
     Real_Array_2D = array of array of real;
     String_Array_2D = array of array of string;
     RealArray3D = array of array of array of real;

     datatype=(real_type,string_type);
     DataFile= record    {this defines a database for the program}
         tvar,  {total variables}
         nvar,  {numeric variables}
         svar,  {string variables}
         nobs,  {number of observations}
         clabellen,  {max length of case labels}
         vlabellen,  {max length of variable labels}
         nvarplace:  {max number of digits to the left of the decimal point, when read}
             integer;
         olabel,     {obs labels, in sdata[*,0]}
         vlabel,     {all variable labels}
         svlabel,    {string variable labels}
         nvlabel:    {numeric variable labels}
             String_Vector;
         vtype:array of datatype;  {type (numeric or string for each varables}
         sdata: String_Array_2D;  {array of string variable values}
                                           {default case labels stored in sdata[0,*]}
         ndata: Real_Array_2D;    {array of numeric variable values}
         svindex,     {svindex[i]=the original variable number of sdata variable i}
         nvindex,     {nvindex[i]=the original variable number of ndata variable i}
         vindex:      {vindex[i]:=the index in sdata(0,-) or ndata(+) of original var i}
           Int_Vector; {array of integer;}
      end;

      DataFile_Int= record    {this defines a database for the program}
         tvar,  {total variables}
         nvar,  {numeric variables}
         svar,  {string variables}
         nobs,  {number of observations}
         clabellen,  {max length of case labels}
         vlabellen,  {max length of variable labels}
         nvarplace:  {max number of digits to the left of the decimal point, when read}
             integer;
         olabel,     {obs labels, in sdata[*,0]}
         vlabel,     {all variable labels}
         svlabel,    {string variable labels}
         nvlabel:    {numeric variable labels}
             String_Vector;
         vtype: array of datatype;  {type (numeric or string for each varables}
         sdata: String_Array_2D;  {array of string variable values}
                                           {default case labels stored in sdata[0,*]}
         ndata: Int_Array_2D;    {array of numeric variable values}
         svindex,     {svindex[i]=the original variable number of sdata variable i}
         nvindex,     {nvindex[i]=the original variable number of ndata variable i}
         vindex:      {vindex[i]:=the index in sdata(0,-) or ndata(+) of original var i}
           Int_Vector; {array of integer;}
      end;
     {end of database definition}

       tokentype=string[tokenlen];
       pathtype=string;

       label_string= string[label_width];
       label_array = array[1..1] of label_string;
       label_pntr  = ^label_array;

       matrix_type = array[1..1] of matrix_elmt;
       matrix_pntr = ^matrix_type;

       TimeDateType=string[8];

{***************************************************************************}

{STRING FUNCTIONS}
function uc(c: ansichar): ansichar;  {convert a char to upper case}

function ucase(strng: ansistring): ansistring;       {convert a string to upper case}

function squash(strng: ansistring): ansistring;     {remove all blanks from a string}

function trim(strng: ansistring): ansistring;    {remove leading and trailing blanks}

function trimR(strng: ansistring): ansistring;   {remove trailing blanks}

function lowstring(xmin: real): tokentype;

function highstring(xmax: real): tokentype;

function padright(strng: ansistring; totallength: integer): ansistring; {pad a string on right}

function truncstring(strng: ansistring; len: integer): ansistring;  {truncate to len}

function intstrfn(i: integer): ansistring;

{INPUT-OUTPUT FUNCTIONS}

function validname(var name: ansistring): boolean;

function file_prefix(name: ansistring): ansistring;   {return name w/o extension}

function file_extension(name: ansistring): ansistring;      {return .ext of file}

function exists(var fvar:text; name: ansistring): boolean;    {does file exist}

function readkey: ansichar;

function readbool(message: ansistring; dflt: ansichar): boolean;
                         {issue message as a prompt and return true or false}

procedure readfile(message: ansistring; var fvar: text; var dflt: ansistring);
    {prompt for an input file name, check its existence and open for reading}

procedure writefile(message: ansistring; var fvar: text; var name: ansistring);
  {prompt for a file name, if it exists, prompt to erase, and open for write}

function readint(message: ansistring; min,max: longint; dflt: ansistring):
  longint;                {return an integer between min and max (inclusive)}

function readreal(message: ansistring; min,max: real; dflt: ansistring): real;
                                         {returns a real; works like readint}

function readchoice(message,choices,dflt: ansistring): ansichar;
         {print message and read one character which must be in choices list}

function readchar(message: ansistring; dflt: ansichar): ansichar;
         {get a nonblank character}

function readline(message,dflt: ansistring): ansistring;

procedure writeplotfile(var plotter: boolean;  var plotsize, plotdest: ansichar;
  var plotfile: text; var plotout: ansistring);

procedure CloseWindow;

procedure ClrScr;

Procedure HaltProgram(msg: ansistring);

{ANTANA FUNCTIONS}

procedure skip_comment(var fvar: text);

function skip_separator(var fvar: text): ansichar;

function gettoken(var fvar: text): tokentype;        {reads an antana token}

function getreal(var fvar: text): real;
                        {reads antana standard reals with embedded comments}

function getinteger(var fvar: text): integer;
                     {reads antana standard integers with embedded comments}

{function getlongint(var fvar: text): longint;}
                     {reads antana standard integers with embedded comments}

function ndx (i,j,cols: integer): integer;
                         {1,1 origin indexing of sequentially stored matrix}

function t_ndx(i,j,size: integer): integer;  {returns triangular matrix indx}

procedure alloc_matrix(var pntr; elmts: longint; size: word);
                                   {generalized memory allocation procedure}

procedure read_antana_matrix(msg: ansistring; var dskin: text;
  var filein: ansistring;  var nrow,ncol: integer;
  var table: matrix_pntr; real_input: boolean);

procedure read_labels(elmt: integer; var maxlen: integer;
  var labl: label_pntr; msg: ansistring; var name: ansistring;
  ext: ansistring; leftjust: boolean);
  {return longest label length}

procedure print_triangle(var mat: matrix_pntr; size: integer ;
  var labl: label_pntr; llen: integer; path,t1,t2,t3: ansistring);

{MATHEMATICAL FUNCTIONS}

function log(val: real): real; { compute log10(val) }

function place(val: real): integer;                { number of digits in val }

function power(base,expon: real):extended;   {raise positive base the expon power}

function normal: extended;           {return a normally distributed random number}

procedure setrandom;

function ceil(x: real): longint;

function floor(x: real): longint;

function maximum(a,b: real): real;

function minimum(a,b: real): real;

function i_maximum(a,b: longint): longint;

function i_minimum(a,b: longint): longint;

function roundceil(r: real): real;

function arcsin(val: real): real;

function twotothe(power: integer): integer;  {use only for integers}

function r_binomial(n,k: longint; p: real): extended;

function lnfactorial(k: longint): extended;

function IsInteger(rv: real; var error: boolean): integer;

function fisher_prob(a,b,c,d: integer; var other_tail: extended): extended;

function chisq_tail_prob(x2:real; dof:integer): real;

function G_compute(var table: Int_Array_2D; var rsum,csum: Int_Vector;
  total: integer): real;

function G_goodfit(var obs: Int_Vector; var expect: Real_Vector;
  vals: integer): real;
  
function williams(var rsum,csum: Int_Vector; total: integer): real;

function williams_goodfit(vals: integer; total: real): real;

{MISCELLANEOUS FUNCTIONS}

procedure sort_real_vector(key: real_vector);  { SORT a column of numbers}

Procedure sort_int_vector(key: int_vector);

procedure boxlinec(len: integer; line: ansistring);

procedure boxline(len: integer; line: ansistring);

procedure boxtop(len: integer; top: boolean);

procedure copyright(name,version,years,description: ansistring);

procedure displayheading(var fvar: text);
  {print # comments in 1st 5 lines of an ADF file}

{TIME AND DATE FUNCTIONS}

function DateString: ansistring;

function ReverseDate(dte: ansistring): ansistring;

function TimeString: ansistring;

function timesec: real;

function readdate(msg: ansistring; dflt: timedatetype): timedatetype;

{CSV Procedures}

function gettoken_csv(var ln: ansistring; var p1,error: integer; var tokentype: datatype): ansistring;
{reads a token from the line}
{note in string type embedded blanks are allowwed but may be incompatible with other programs}

function get_csv(var f: datafile; prompt: ansistring; var filename: ansistring; quiet: boolean): boolean;

function ReadFile_CSV_Int(var dfi: DataFile_Int; prompt: ansistring; var filename: ansistring): boolean;

procedure readfile_csv(var df: datafile; prompt: ansistring; var filename: ansistring);

procedure readfile_csv_quiet(var df: datafile; prompt: ansistring; var filename: ansistring);

procedure datafiletoscreen_csv(var f: datafile);

procedure datafiletoprinter_csv(var f:datafile; prompt: ansistring);

function csvstring(s: ansistring): ansistring;

function csvreal(v: extended; w,d: integer): ansistring;

procedure datafiletofile_csv(var f: datafile; prompt: ansistring);


{MATRIX UTILITIES}

function vector_median(var vect: Real_Vector; elmt:integer): real;

procedure Clear_2D_Int_Array(var table: Int_Array_2D);

procedure Fill_Work(var itable: Int_Array_2D; var rtable: Real_Array_2D);

function Int_Vector_Is_0(var vector: Int_Vector): boolean;

procedure Multiply_2D_Real_Array(var table: Real_Array_2D; factor: real);

procedure Multiply_Real_Row_Or_Col(var table: Real_Array_2D; rowcol: integer;
  rowwise: boolean; factor: real);

procedure Multiply_Real_Vector(var vector: Real_Vector; factor: real);

function Real_Row_Max(var table: Real_Array_2D; rowno: integer): real;

function Real_Column_Max(var table: Real_Array_2D; colno: integer): real;

function Real_Vector_Max(var vector: Real_Vector): real;

function Real_Vector_Min(var vector: Real_Vector): real;

function Int_Vector_Min(var vector: Int_Vector; elmts: integer): longint;

function String_Vector_Maxlen(var s: String_Vector): integer;

procedure Margins_2D_Int_Array(var table: Int_Array_2D; var rowsum,colsum: Int_Vector;
  var total: integer);

procedure Margins_2D_Real_Array(var table: Real_Array_2D; var rowsum,colsum: Real_Vector;
  var total: real);

procedure pause(var fle: ansistring);

{*************************************************************************}

implementation

uses SysUtils;

{STRING FUNCTIONS}

function uc(c: ansichar): ansichar;
var n: byte;
begin
  {n:=ord(c);
  if (n>+97) and (n<=122) then c:=chr(n-32);
  uc:=c;}
  uc:=upcase(c);
end;

function ucase(strng: ansistring): ansistring;
{change all lower case letters in strng to upper case}
var i,j,n: integer;
begin
  {j:=length(strng);
  for i:=1 to j do begin
    n:=ord(strng[i]);
    if (n>=97) and (n<=122) then strng[i]:=chr(n-32);
  end;
  ucase:=strng; }
  ucase:=UpperCase(strng);
end;

function squash(strng: ansistring): ansistring;
{eliminate all blanks from strng}
var j: integer;
begin
  j:=length(strng);
  while j>=1 do begin
    if strng[j]=' ' then
      strng:=copy(strng,1,j-1)+copy(strng,j+1,length(strng)-j);
    j:=j-1;
  end;
  squash:=strng;
end;

function trim(strng: ansistring): ansistring;
{eliminate leading and trailing blanks from strng}
var i,j: integer;
begin
  i:=length(strng);
  while (i>0) and (strng[i]=' ') do i:=i-1;
  strng:=copy(strng,1,i);
  j:=1;
  while (J<I) and (strng[j]=' ') do j:=j+1;
  trim:=copy(strng,j,i-j+1);
end;

function trimR(strng: ansistring): ansistring;
{eliminate trailing blanks from strng}
var i: integer;
begin
  i:=length(strng);
  while (i>0) and (strng[i]=' ') do i:=i-1;
  trimr:=copy(strng,1,i);
end;

function lowstring(xmin: real): tokentype;
var strng: tokentype;
begin
  str(floor(xmin):16,strng);
  strng:=trim(strng);
  lowstring:=strng;
end;

function highstring(xmax: real): tokentype;
var strng: tokentype;
begin
  str(ceil(xmax):16,strng);
  strng:=trim(strng);
  highstring:=strng;
end;

function padright(strng: ansistring; totallength: integer): ansistring; {pad a string on right}
var i: integer;
begin
  for i:=length(strng)+1 to totallength do strng:=strng+' ';
  padright:=strng;
end;

function truncstring(strng: ansistring; len: integer): ansistring;  {truncate to len}
begin
  if length(strng)<=len then truncstring:=strng else truncstring:=copy(strng,1,len);
end;

function intstrfn(i: integer): ansistring;
var rslt: ansistring;
begin
  str(i,rslt);
  IntStrFn:=rslt;
end;

{INPUT-OUTPUT FUNCTIONS}

function validname(var name: ansistring): boolean;
{check to see if name is a valid pathname, return true or false}
{this procedure is not perfect}
    {i.e. it will let some invalid names by, but most will be caught}
const illegal='/*?"<>|';
var ok,logicaldevice: boolean; i,len: integer; c: ansichar;
begin
  len:=length(name);
  name:=ucase(name);

  logicaldevice:=
        ((length(name)=3) and (pos(name,'PRN|NUL')<>0)) or
        ((length(name)=4) and (pos(name,'COM1|COM2|LPT1|LPT2|LPT3')<>0));
  i:=1;
  ok:=not logicaldevice;
  while i<=len do begin
    c:=name[i];
    if (pos(c,illegal)>0) or (ord(c)<32) or (ord(c)>=127) then ok:=false;
    inc(i);
  end;
  validname:=ok and (len>0);
end;

function file_prefix(name: ansistring): ansistring; {return name w/o extension}
var len,pos: integer;
begin
  len:=length(name);
  pos:=len;
  while (pos>1) and (name[pos]<>'.') and (pos>(len-3))do dec(pos);
  if (pos>0) and (name[pos]='.') then file_prefix:=copy(name,1,pos-1)
  else file_prefix:=name;
end;

function file_extension(name: ansistring): ansistring; {return extension including '.'}
var len,pos: integer;
begin
  len:=length(name);
  pos:=len;
  while (pos>1) and (name[pos]<>'.') and (pos>(len-3))do dec(pos);
  if name[pos]='.' then file_extension:=copy(name,pos+1,len-pos+1)
  else file_extension:='';
end;

function exists(var fvar:text; name: ansistring): boolean;
var ok: boolean;
begin
  AssignFile(fvar,name);
  {$I-}  reset(fvar);  {$I+}
  ok:=(ioresult=0);
  if ok then close(fvar);
  exists:=ok;
end;

function readkey: ansichar;
{reads the rest of the line and extracts the first character}
{if the line is empty returns <CR>}
{replacement for DOS readkey which does not exist in Delphi}
var c: ansichar; s: ansistring;
begin
  readln(s);
  if length(s)=0 then c:=chr(13)
  else c:=s[1];
  readkey:=c;
end;

function readbool(message: ansistring; dflt: ansichar): boolean;
{issue message as a prompt and return true or false}
  {the procedure reprompts until a valid response is received }
  {if dflt='T' then a <CR> will result in true being returned }
  {if dflt='F' then a <CR> will result in false being returned}
  {otherwise the program requires a definitive response       }
  {y,ye,yes,n,no,ok (any case) are accepted as valid responses}
var ok: boolean;  reply,c: ansichar;
begin
  ok:=true;  c:=' ';  readbool:=true;
  repeat
    {dflt 'T' true; 'F' false; otherwise no default}
    if ok then begin
      if dflt='T' then c:='Y' else if dflt='F' then c:='N' else c:=' ';
    end;
    if c=' ' then write(message,' ? ')
    else write(message,' {',c,'} ? ');
    reply:=readkey;
    if (reply=#13) or (reply=#10) then reply:=c;
    reply:=upcase(reply);
    ok:=true;
    case reply of
      'Y': readbool:=true;
      'N': readbool:=false;
      else begin
             ok:=false;
             {write(#7);}
           end;
    end;
  until ok;
  {writeln(reply);}
end;

procedure readfile(message: ansistring; var fvar: text; var dflt: ansistring);
{prompt for an input file name, check its existence and open for reading}
var  ok: boolean;  name: ansistring;
begin
  {if length(dflt)=0 then dflt:='*';}
  ok:=false;
  if ((pos('.',name)=4) and (pos(copy(name,1,3),'CON|PRN|NUL')<>0)) or
     ((pos('.',name)=5) and (pos(copy(name,1,4),'COM1|COM2|LPT1|LPT2|LPT3')<>0))
    then name:=copy(name,1,pos('.',name)-1);

  while not OK do begin
    if dflt='' then write(message,' ? ') else write(message,' {',dflt,'} ? ');
    readln(name);
    name:=trim(name);
    {if name='?' then directory;}
    if (length(name)=0) and (length(dflt)>0) and (dflt[1]<>'.') then
      name:=dflt;
    ok:=validname(name);
    if ok and (name<>'CON') and ((length(dflt)>0) and (dflt[1]='.')) and (pos('.',name)=0) then begin
      name:=name+dflt;  ok:=validname(name); {recheck w/ added ext}
    end;
    if ok then begin
      if name='CON' then assign(fvar,'')
      else assign(fvar,name);
      {$I-}  reset(fvar);  {$I+}
      ok:=(ioresult=0);
      {ok:=ok and (pos(name,'CON|')<>0);}
      if not ok then writeln('File Not Found');
    end;
  end;
  dflt:=name;
end;

procedure writefile(message: ansistring; var fvar: text; var name: ansistring);
{prompt for a file name, if it exists, prompt to erase, and open for write}
var
  ok: boolean; dflt: ansistring;
begin
  {ok:=false;}

  if ((pos('.',name)=4) and (pos(copy(name,1,3),'CON|PRN|NUL')<>0)) or
     ((pos('.',name)=5) and (pos(copy(name,1,4),'COM1|COM2|LPT1|LPT2|LPT3')<>0))
    then name:=copy(name,1,pos('.',name)-1);

  dflt:=name;

  {if length(dflt)=0 then dflt:='*';}
  repeat
    if dflt='' then write(message,' ? ')
    else write(message,' {',dflt,'} ? ');
    readln(name);
    name:=trim(name);
    { if name='?' then directory; }
    if (length(name)=0) and (length(dflt)>0) and (dflt[1]<>'.') then
      name:=dflt;
    ok:=validname(name);
    if ok and (name<>'CON') and ((length(dflt)>0) and (dflt[1]='.'))  and (pos('.',name)=0) then begin
      name:=name+dflt;  ok:=validname(name); {recheck w/ added ext}
    end;
    if ok then begin
      if name='CON' then AssignFile(fvar,'')
      else if exists(fvar,name) then begin
        ok:=readbool('File Exists; OK to Erase It','F');
        if ok then begin
          erase(fvar);  writeln(name,' Erased');
        end;
      end
      else ok:=true;
      if ok then begin
        rewrite(fvar);
        ok:=(ioresult=0);
        if not ok then writeln(' I/O Error - Disk May Be Full');
      end;
    end;
  until ok;
end;

function readint(message: ansistring; min,max: longint; dflt: ansistring):
  longint;
{return an integer between min and max (inclusive)}
{if the dflt is specified AS A STRING, it will be returned if <CR> is entered}
var ok: boolean; reply: ansistring; value: longint; error: integer;
begin
  ok:=false; readint:=0;
  repeat
    if dflt<>'' then write(message,' {',dflt,'} ? ')
    else write(message,' ? ');
    readln(reply);
    reply:=trim(reply);
    if (reply<>'') or (dflt<>'') then begin
      if length(reply)=0 then reply:=dflt;
      val(reply,value,error);
      ok:=(error=0);
      if not ok then writeln('Invalid Character in Integer: "',
        reply[error],'"')
      else begin
        ok:=(value>=min) and (value<=max);
        if not ok then
          writeln('Value Out of Range (',min,',',max,')')
        else readint:=value;
      end;
    end;
  until ok;
end;

function readreal(message: ansistring; min,max: real; dflt: ansistring): real;
{returns a real; works like readint}
var ok: boolean; reply: ansistring; value: real; error: integer;
begin
  ok:=false;  readreal:=0.0;
  repeat
    if dflt<>'' then write(message,' {',dflt,'} ? ')
    else write(message,' ? ');
    readln(reply);
    reply:=trim(reply);
    if (reply<>'') and (reply[1]='.') then reply:='0'+reply
    else if (reply<>'') and (reply[1]='-') and (reply[2]='.') then
      reply:='-0.'+copy(reply,3,length(reply)-2);
    if (reply<>'') or (dflt<>'') then begin
      if length(reply)=0 then reply:=dflt;
      val(reply,value,error);
      ok:=(error=0);
      if not ok then writeln('Invalid Character in Real: "',reply[error],'"')
      else begin
        ok:=(value>=min) and (value<=max);
        if not ok then
          writeln('Value Out of Range (',min,',',max,')')
        else readreal:=value;
      end;
    end;
  until ok;
end;

function readchoice(message,choices,dflt: ansistring): ansichar;
{print message and read one character which must be in choices list}
{alpha characters in choices must be upper case, default char=dflt}
var ok: boolean;  reply: ansichar;
begin
  {ok:=true;}
  repeat
    {if ok then} begin
      if dflt='' then write(message,' ? ')
      else write(message,' {',dflt,'} ? ');
    end;
    reply:=readkey;
    if (reply=#13) or (reply=#10) then begin
      if dflt<>'' then reply:=dflt[1] else reply:=chr(255);
    end;
    reply:=upcase(reply);
    ok:=false;
    if pos(reply,choices)>0 then ok:=true else write(#7);
  until ok;
  {writeln(reply);}
  readchoice:=reply;
end;

function readchar(message: ansistring; dflt: ansichar): ansichar;
{get a blank character, if 1st char of message=chr(0), char may be blank}
var ok,allowblank: boolean;  reply,c: ansichar; minchar: byte;
begin
  {ok:=true;}
  allowblank:= ord(message[1])=0;
  if allowblank then begin
    minchar:=32;
    message:=copy(message,2,length(message)-1);
  end
  else minchar:=33;
  repeat
    {if dflt=' ' or nul then no default}
    {if ok then} begin
      if dflt=chr(0) then c:=' ' else c:=dflt;
      if c=' ' then write(message,' ? ')
      else write(message,' {',c,'} ? ');
    end;
    reply:=readkey;

    if (reply=#13) or (reply=#10) then reply:=c;
    ok:=true;
    if ord(reply)<minchar then begin
      if reply=chr(0) then reply:=readkey; {flush scan code}
      ok:=false;
      write(#7);
    end;
  until ok;
  {writeln(reply);}
  readchar:=reply;
end;

function readline(message,dflt: ansistring): ansistring;
var reply: ansistring;
begin
  if length(dflt)>0 then
    write(message,'{',dflt,'} ? ')
  else write(message,' ? ');
  readln(reply);
  if reply='' then reply:=dflt;
  readline:=reply;
end;

procedure writeplotfile(var plotter: boolean; var plotsize, plotdest: ansichar;
  var plotfile: text; var plotout: ansistring);
var dflt: string[1];  {5.1}
begin
  plotter:=readbool('Plot to HPGL Plotter or HP LaserJet Printer','F');
  if plotter then begin
    plotout:='';
    writefile('HPGL Output to (PRN, COM1, COM2 or <filename>)',plotfile,
      plotout);
    if plotout='PRN' then dflt:='L' else dflt:='A';
    plotsize:=readchoice(
      ' HP [L]aserJet 3/4 or Plotter Paper Size: [A] 8.5x11 [B] 11x17',
      'ABL',dflt);
    plotout:=trim(plotout);
    plotout:=ucase(plotout);
    if (plotout='COM1') or (plotout='COM2') or (plotout='PRN') or
      (plotout='LPT1') or (plotout='LPT2') then plotdest:='D'
      else plotdest:='F';
    if plotsize='L' then begin
      plotsize:='A';
      {plotdest:=chr(ord(plotdest)+32);}  {sets to d or f instead of D or F}
      if plotdest='D' then plotdest:='d' else if plotdest='F' then plotdest:='f';
    end;
    if (plotdest='D') and (plotsize<>'L') then begin
      writeln('  Plotter Must Be Turned on and Connected to ',plotout);
      if plotout[1]='C' then
        writeln('  Port Must Be Initialized with MODE ',plotout,':baud,,,,P');
      if readbool('  Abort Now','F') then
        HaltProgram('Plotter Not Initialized');
    end;
  end else plotsize:='A';
end;

procedure CloseWindow;
begin
  repeat
  until readbool('OK to Close Program Window','T');
end;

Procedure ClrScr;
var i: byte;
begin
  for i:=1 to 25 do writeln;
end;

procedure HaltProgram(msg: ansistring);
begin
  writeln;
  writeln('<<<<<<<<<<<<>>>>>>>>>>>>');
  writeln(msg,' - Program Aborting.');
  CloseWindow;
  Halt;
end;


{ANTANA FUNCTIONS}

procedure skip_comment(var fvar: text);
var c: ansichar;
begin
  {skip over # to end of line or text between #'s, e.g. #site 12 #}
  c:=^J;
  repeat
    if not seekeoln(fvar) then read(fvar,c);
  until eoln(fvar) or (c='#') or (c=^J);
end;

function skip_separator(var fvar: text): ansichar;
var c: ansichar;
begin
  {skip over blanks, tabs, commas, lf's, comments}
  c:=^Z;
  repeat
    if not seekeof(fvar) then begin
      read(fvar,c);
      if c='#' then begin skip_comment(fvar); c:=',' end;
    end;
  until eof(fvar) or (pos(c,',#')=0);
  skip_separator:=c;
end;

function gettoken(var fvar: text): tokentype;
{reads an antana token}
var reply: string[32]; len: integer; c: ansichar;
begin
  c:=skip_separator(fvar);
  len:=0; reply:='';
  repeat
    len:=succ(len);
    reply:=reply+c;
    if eoln(fvar) {or eof(fvar)} then c:=' ' else read(fvar,c);
  until pos(c,' ,#'^I^J)>0;
  if c='#' then skip_comment(fvar);
  {reply[0]:=chr(len);}
  gettoken:=leftstr(reply,len);
  {gettoken:=reply;}
end;

function getreal(var fvar: text): real;
{reads antana standard reals with embedded comments}
var reply: string[32]; value: real; error: integer;
begin
  reply:=gettoken(fvar);
  getreal:=missing;
  if ord(reply[1])=26 then writeln('Unexpected EOF')
  else if pos('?',reply)=0 then begin
    val(reply,value,error);
    if error<>0 then begin
      writeln('Invalid Character in Number: "',reply[error],'"',
        '  (Decimal ',ord(reply[error]),') - Input: ',reply);
      if readbool('Quit Now','T') then begin
        close(fvar); HaltProgram('Real Input Error');
      end;
    end
    else getreal:=value;
  end;
end;

function getinteger(var fvar: text): integer;
{reads antana standard long integer with embedded comments}
var rvalue: real;
begin
  getinteger:=-maxint;
  rvalue:=getreal(fvar);
  if (rvalue>high(integer)) or (rvalue<low(integer)) or ((rvalue-int(rvalue))<>0.0) then
  begin
    writeln('Not an Integer or Integer Out of Range: ',rvalue:16:4);
    if readbool('Quit Now','T') then begin
      close(fvar); HaltProgram('Integer Input Error');
    end;
  end
  else getinteger:=round(rvalue);
end;

{function getlongint(var fvar: text): longint;}
{KWK 12/03 Eliminated because integer is now same as longint}
{reads antana standard long integer with embedded comments}
{var rvalue: real;
begin
  getlongint:=-maxlongint;
  rvalue:=getreal(fvar);
  if (abs(rvalue)>maxlongint) or ((rvalue-int(rvalue))<>0.0) then
  begin
    writeln('Not an Integer or Integer Out of Range: ',rvalue:16:4);
    if readbool('Quit Now','T') then begin
      close(fvar); HaltProgram('Integer Input Error');
    end;
  end
  else getlongint:=round(rvalue);
end;}

function ndx(i,j,cols: integer): integer;  {1,1 origin indexing}
begin
  ndx := (i-1)*cols+j;
end;

function t_ndx(i,j,size: integer): integer;
var rslt,t: integer;  {no diagonal}
begin
  i:=pred(i);  j:=pred(j);
  if i>j then begin
    t:=i;  i:=j;  j:=t;
  end;
  rslt:=j-i+((2*size-i-1)*i) div 2; {assumes i<j}
  t_ndx:=rslt;
end;

procedure alloc_matrix(var pntr; elmts: longint; size: word);
var bytes: longint;
begin
  bytes:=elmts*size;
  try
    getmem(pointer(pntr),bytes);
  except
    on EOutOfMemory do begin
      writeln('Insufficient Memory to Allocate Matrix; Program Halts');
      CloseWindow;
      halt;
    end;
  end;
end;
{$R-}
procedure read_antana_matrix(msg: ansistring; var dskin: text;
  var filein: ansistring;  var nrow,ncol: integer;
  var table: matrix_pntr; real_input: boolean);
var i,j: integer;
begin
  readfile(msg,dskin,filein);
  if filein='CON' then begin
    nrow:=readint('  Number of Rows (Observations)',1,maxint,'');
    ncol:=readint('  Number of Columns (Variables)',1,maxint,'');
    writeln('  Enter ',nrow*ncol,' Values Followed by <Enter>');
  end
  else begin
    displayheading(dskin);
    nrow:=getinteger(dskin);  ncol:=getinteger(dskin);
    writeln('  Reading ',nrow,' Rows and ',ncol,' Columns');
  end;
  alloc_matrix(table,nrow*ncol,sizeof(matrix_elmt));
  if filein='CON' then write('? ');
  for i:=1 to nrow do
    for j:=1 to ncol do begin
      if (filein='CON') and (eoln(dskin)) then write('? ');
      if real_input then table^[ndx(i,j,ncol)]:=getreal(dskin)
      else table^[ndx(i,j,ncol)]:=getinteger(dskin);
    end;
  close(dskin);
end;

procedure read_labels(elmt: integer; var maxlen: integer;
  var labl: label_pntr; msg: ansistring;  var name: ansistring;
  ext: ansistring; leftjust: boolean);
var uno,ll: integer;  strng: label_string;  dskin: text;
begin
  if readbool('Read '+msg+' Label File','F') then begin
    name:=file_prefix(name)+'.'+ext;
    readfile('  '+msg+' Label File Name',dskin,name);
    for uno:=1 to elmt do begin
      readln(dskin,strng);
      labl^[uno]:=trim(strng);
      if length(labl^[uno])>maxlen then maxlen:=length(labl^[uno]);
    end;
    close(dskin);
    if leftjust then for uno:=1 to elmt do
      for ll:=length(labl^[uno])+1 to maxlen do labl^[uno]:=labl^[uno]+' ';
  end else
    for uno:=1 to elmt do if leftjust then str(uno:maxlen,labl^[uno])
      else str(uno,labl^[uno]);
end;

procedure print_triangle(var mat: matrix_pntr; size: integer ;
  var labl: label_pntr; llen: integer; path,t1,t2,t3: ansistring);
var  i,hlen,linesize,itemwid,dec,maxplace,elmts: integer;
     triangle: boolean;
     numstr: string[3];
     fileout: ansistring;
     dskout: text;
     scale: real;

procedure print_tri_matrix;
var i,j,ss,across,subset,llim,ulim: integer;
begin
  across:=(linesize-llen-1) div itemwid;
  subset:=(size+across-1) div across-1;
  writeln(dskout);
  if length(t1)>0 then writeln(dskout,t1);
  if length(t2)>0 then writeln(dskout,t2);
  if length(t3)>0 then writeln(dskout,t3);
  for ss:=0 to subset do begin
    if triangle then llim:=ss*across+2 else llim:=1;
    for i:=llim to size do begin
      write(dskout,labl^[i]:llen,':');
      if triangle then ulim:=i-1 else ulim:=size;
      for j:=ss*across+1 to i_minimum(succ(ss)*across,ulim) do begin
        if i<j then write(dskout,mat^[t_ndx(i,j,size)]*scale:itemwid:dec)
        else if i=j then write(dskout,' ':itemwid)
        else write(dskout,mat^[t_ndx(j,i,size)]*scale:itemwid:dec);
      end;
      writeln(dskout);
    end;
    write(dskout,' ':llen+1);
    if triangle then ulim:=pred(size) else ulim:=size;
    for j:=ss*across+1 to i_minimum(succ(ss)*across,ulim) do
      write(dskout,' ',copy(labl^[j],1,itemwid-1):itemwid-1);
    writeln(dskout);
    writeln(dskout);
  end;
end;

begin
  elmts:=size*pred(size) div 2;
  maxplace:=0;
  for i:=1 to elmts do begin
    {write(mat^[i]:8);}
    maxplace:=i_maximum(maxplace,place(mat^[i]));
  end;

  fileout:=file_prefix(path)+'.TXT';
  writefile('Output File Name',dskout,fileout);

  linesize:=80;
  if maxplace>2 then dec:=0 else dec:=2;
  scale:=1.0;
  hlen:=llen;
  triangle:=true;

  if readbool('Set Output Options','F') then begin
    triangle:=readchoice('  Print Matrix in [T]riangular or [S]quare Form',
      'TS','T')='T';
    str(linesize,numstr);
    linesize:=readint('  Output Line Width (in characters)',40,maxint,numstr);
    str(hlen,numstr);
    hlen:=readint('  Print Width of Observation Labels',3,label_width,numstr);
    scale:=readreal('  Scale Factor (Multiply Values by)',0.001,1000,'1');
    if scale>=10.0 then dec:=i_minimum(0,dec-place(scale)+1);
    str(dec,numstr);
    dec:=readint('  Number of Decimal Places to Print',0,16,numstr);
  end;

  itemwid:=i_maximum(succ(hlen),dec+maxplace+2);

  print_tri_matrix;
  close(dskout);
end;

{MATHEMATICAL FUNCTIONS}

function log(val: real): real; { compute log10(val) }
begin
  log:=ln(val)/ln(10);
end;

function place(val: real): integer;  {number of digits left of decimal in val}
var neg: integer;
begin
  if val<0.0 then begin
    neg:=1; val:=abs(val);
  end else neg:=0;
  if val<1.0 then val:=1.0;
  place:=succ(trunc(log(val+0.0000001)))+neg;
end;

function power(base,expon: real):extended;  {raise base the expon power}
{works only for positive base}
begin
  if base=0.0 then power:=0.0
  else power:=exp(expon*ln(base));
end;

function normal: extended;
var recalc: boolean;
    r1, r2: real;
{const recalc: boolean=true;
      r1: real=0.0;
      r2: real=0.0;}  {Delphi xe2 Change}

procedure nran;
var x,y,xy,z: real;
begin
  repeat
    x:=2.0*random-1.0;
    y:=2.0*random-1.0;
    xy:=sqr(x)+sqr(y);
  until xy<1.0;
  z:=sqrt(-2.0*ln(xy)/xy);
  r1:=x*z;
  r2:=y*z;
end;

begin
  recalc:=true;
  r1:=0.0; r2:=0.0;
  if recalc then begin
    nran;
    normal:=r1;
  end
  else normal:=r2;
  recalc:=not recalc;
end;

procedure setrandom;
begin
  randseed:=readint('Random Generator Seed (0 to set from clock)',
    -maxlongint,maxlongint,'0');
  if randseed=0 then randomize;
end;

function ceil(x: real): longint;
begin
  
  if (frac(x)=0.0) or (x<0.0) then ceil:=trunc(x)
  else if x>0 then ceil:=trunc(x+1.0)else ceil:=0;
end;

function floor(x: real): longint;
begin

  if (frac(x)=0.0) or (x>0) then floor:=trunc(x)
  else floor:=trunc(x-1.0);
end;

function maximum(a,b: real): real;
begin
  if a>b then maximum:=a else maximum:=b;
end;

function minimum(a,b: real): real;
begin
  if a>b then minimum:=b else minimum:=a;
end;

function i_maximum(a,b: longint): longint;
begin
  if a>b then i_maximum:=a else i_maximum:=b;
end;

function i_minimum(a,b: longint): longint;
begin
  if a>b then i_minimum:=b else i_minimum:=a;
end;

function roundceil(r: real): real;
var t: real; i: integer;
begin
  if r<0.0 then begin
    if r>=-10.0 then roundceil:=ceil(r)
    else begin
      t:=power(10,ceil(log(-r))-1.0);
      i:=-10;
      while (t*i)<=r do i:=i+1;
      roundceil:=t*i;
    end;
  end
  else begin
    if r<=10.0 then roundceil:=ceil(r)
    else begin
      t:=power(10,ceil(log(r))-1.0);
      i:=1;
      while (t*i)<r do i:=i+1;
      roundceil:=t*i;
    end;
  end;
end;

function arcsin(val: real): real;
begin
  if abs(val-1.0)<1e-10 then val:=arctan(1e10) else
  val:=arctan(val/sqrt(1.0-sqr(val)));
  arcsin:=val;
end;

function twotothe(power: integer): integer;
begin
  if (power>=0) and (power<=31) then twotothe:=round(exp(power*ln(2)))
  else twotothe:=0;
end;

function r_binomial(n,k: longint; p: real): extended;
{recursive calculation of binomial probability}
{cumulative probabilities computed w/o extra work; tail probs include k}
var r: longint; ltail, rtail, b: extended; q, expect, obs: real;
begin
  r:=0;
  q:=1.0-p;
  if p>=1.0 then begin
    if k<n then begin
      r_binomial:=0.0;  ltail:=0.0;  rtail:=1.0;
    end
    else begin
      r_binomial:=1.0;  ltail:=1.0;  rtail:=1.0;
    end;
  end else
  begin
    B:=power(q,n);             {B(0)=q^n}
    if B=0.0 then writeln('Cannot be Computed Recursively') {underflow}
    else begin
      ltail:=B;
      while r<k do begin
        inc(r);                  {recursive rule}
        B:=B*p*(n-r+1)/(q*(r));  {B(r)=B(r-1)*p*(n-r+1)/(q*r)}
        ltail:=ltail+B;
      end;
    end;
    rtail:=1.0-ltail+B;
    expect:=p*n;  obs:=k;
    if expect=obs then r_binomial:=1.0 {special case, subject to representation error}
    else if expect>obs then r_binomial:=-ltail else r_binomial:=rtail;
  end;
end;

function lnfactorial(k: longint): extended;
var i: longint; fact: extended;
begin
  if k<0 then fact:=0 else if k<=1754 {10} then begin
    fact:=1.0;
    for i:=2 to k do fact:=fact*i;
    fact:=ln(fact);
  end
  else fact:=(k+0.5)*ln(k)-k+0.918938534+1.0/(12.0*k);
  {ln(k!)=(k+.5)ln(k)-k+.5pi+1/12k Stirling's Formula}
  lnfactorial:=fact;
end;

function IsInteger(rv: real; var error: boolean): integer;
begin
  error:=false;
  if (frac(rv)=0.0) and (int(rv)<=high(integer)) and (int(rv)>=low(integer)) then
    IsInteger:=trunc(rv)
  else begin
    IsInteger:=0;  error:=true;
  end;
end;

function fisher_prob(a,b,c,d: integer; var other_tail: extended): extended;
  function constant_part(a,b,c,d: integer): extended;
  begin
    constant_part:=lnfactorial(a+b)+lnfactorial(c+d)+lnfactorial(b+d)+
      lnfactorial(a+c)-lnfactorial(a+b+c+d);
  end;

  function exact_prob(a,b,c,d: integer; term: extended): extended;
  begin
    exact_prob:=exp(term-(lnfactorial(a)+lnfactorial(b)+lnfactorial(c)+
            lnfactorial(d)));
  end;

var total_prob,constant_term,tail_prob: extended;
    oa,ob,oc,od: integer;

begin
  {total_prob:=0.0;}
  {constant_term:=0.0;}
  oa:=a;  ob:=b;  oc:=c; od:=d;

  constant_term:=constant_part(a,b,c,d);

  {Compute First Tail}
  total_prob:=exact_prob(a,b,c,d,constant_term);
  while (a>0) and (b>0) and (c>0) and (d>0) do begin
    if (a*d-b*c)<0 then begin
      dec(a);  dec(d);  inc(b);  inc(c);
    end
    else begin
      dec(b);  dec(c);  inc(a);  inc(d);
    end;
    total_prob:=total_prob+exact_prob(a,b,c,d,constant_term);
  end;

  {compute other tail}
  a:=oa;  b:=ob;  c:=oc;  d:=od;
  if (a*d-b*c)<0 then begin
    if b<c then begin
      a:=a+b;  c:=c-b;  d:=d+b;  b:=0;
    end
    else begin
      a:=a+c;  b:=b-c;  d:=d+c;  c:=0;
    end;
  end
  else begin
    if a<d then begin
      b:=b+a; c:=c+a;  d:=d-a;  a:=0;
    end
    else begin
      a:=a-d;  b:=b+d;  c:=c+d;  d:=0;
    end;
  end;
  other_tail:=0.0;
  tail_prob:=exact_prob(a,b,c,d,constant_term);
  while ((tail_prob+other_tail)<total_prob) and
    ((a<>oa) or (b<>ob) or (c<>oc) or (d<>od)) do begin
    other_tail:=other_tail+tail_prob;
    if (a*d-b*c)>0 then begin
      dec(a);  dec(d);  inc(b);  inc(c);
    end
    else begin
      dec(b);  dec(c);  inc(a);  inc(d);
    end;
    tail_prob:=exact_prob(a,b,c,d,constant_term);
  end;
  other_tail:=other_tail+total_prob;

  fisher_prob:=total_prob;
end;

function chisq_tail_prob(x2:real; dof:integer): real;
{** chisq procedure declarations ** from Henry Harpending **}
  function gfunc(t:real): real;
  var  g:real; j: integer;
  begin
    if (t=0) then g := 1.0
    else if ( abs(1.0-t) > 0.1) then
      g := ( (1.0 - t*t) + 2.0 * t * ln(t))/ ( (1.0-t) * (1.0-t) )
    else begin
      g := 0.0;
      for j := 1 to 5 do  {note KK added Abs to power below}
        g := g + 2.0 * (power(abs(1.0-t) , j)) / ( (1.0+j)*(2.0+j) );
    end;
    gfunc := g;
  end;

  function fnl(a:real):real;
  begin
    fnl:=(1.0/(12.0*a))*(1.0-(1.0/a)*(1.0/30.0-(1.0/a)*
      (1.0/105.0-1.0/(140.0*a))));
  end;

  function sign(x: real): real;
  begin
    if(x<0.0) then sign := -1.0 else sign := 1.0;
  end;

  function normalarea(z:real):real;
  var  t:real;
  begin
    t := 0.806*abs(z)*(1.0 - 0.018*abs(z));
    normalarea := 0.5 + 0.5 * sign(z) * sqrt(1.0 - exp(- t*t));
  end;

var  {j, result: integer;}
     t2,z,z1,{z2,}c,g,d,a,d3,q1{,d1}: real;
begin
  if ( dof<11 ) then begin
    if x2=0.0 then q1:=1.0 else
    if x2>1000 then q1:=0.0
    else begin
      z := x2/2.0;   {z2 := z*z;}    c := 1.0;    g := 1.0;
      d := dof/2.0;    a := d;     d3 := d+2.0;
      repeat
        a := a+1.0;   c := c*z/a;   g := g+c;
      until ( (c/g) < 0.5e-6);
      g := g * exp (d * ln(z) - d3 * fnl(d3*d3) - (d3 - 0.5) * ln(d3) +d3 -z) * (d+1.0);
      q1 := 1.0 -g/sqrt(2.0*3.14159);
    end;
  end
  else
  begin {equivalent normal deviate}
    t2 := (dof-1.0)/x2;
    d3 := (x2 -dof) + (0.6666667 -0.08/dof);
    z1 := d3 * sqrt((1.0+gfunc(t2)) / (2.0*x2));
    q1 := 1.0-normalarea(z1);
  end;
  chisq_tail_prob:=q1;
end;


function G_compute(var table: Int_Array_2D; var rsum,csum: Int_Vector;
  total: integer): real;
var i,j,rows,cols: integer; gsum: real;
{assumes marginals computed}
begin
  rows:=high(table);  cols:=high(table[0]);
  gsum:=0.0;
  for i:=1 to rows do for j:=1 to cols do
    if table[i,j]>0 then gsum:=gsum+table[i,j]*ln(table[i,j]);
  for i:=1 to rows do gsum:=gsum-rsum[i]*ln(rsum[i]);
  for j:=1 to cols do gsum:=gsum-csum[j]*ln(csum[j]);
  gsum:=gsum+total*ln(total);
  G_compute:=2.0*gsum;
end;

function G_goodfit(var obs: Int_Vector; var expect: Real_Vector;
  vals: integer): real;
var i: integer; gsum: real;
{assumes marginals computed}
begin
  gsum:=0.0;
  for i:=1 to vals do
    if obs[i]>0 then gsum:=gsum+obs[i]*ln(obs[i]/expect[i]);
  G_goodfit:=2.0*gsum;
end;

function williams(var rsum,csum: Int_Vector; total: integer): real;
var i,j,rows,cols: integer; rf,cf,q: real;
{assumes marginals computed}
begin
  rows:=high(rsum);  cols:=high(csum);
  rf:=0.0;
  for i:=1 to rows do rf:=rf+1.0/rsum[i];
  rf:=total*rf-1.0;
  cf:=0.0;
  for j:=1 to cols do cf:=cf+1.0/csum[j];
  cf:=total*cf-1.0;
  q:=1+rf*cf/(6.0*total*(rows-1)*(cols-1));
  williams:=q;
end;

function williams_goodfit(vals: integer; total: real): real;
begin
  williams_goodfit:=1.0+(vals+1)/(6*total);
end;


{MISCELLANEOUS FUNCTIONS}

procedure sort_real_vector(key: real_vector);  { SORT a column of numbers}
  var n,part,h,i,j,s: integer;  k,t: real;  stop: boolean;
begin  { quicksort see knuth Vol 3 }
  n:=high(key);
  part:=1;
  while (twotothe(part+1)-1)<(n div 3) do part:=part+1;
  for s:=part downto 1 do begin
    h:=twotothe(s-1);
    for j:=h+1 to n do begin
      stop:=false;  i:=j-h;  k:=key[j];
      while (i>0) and not stop do
        if k<key[i] then begin
          t:=key[i]; key[i+h]:=t;
          {key[i+h]:=key[i];} i:=i-h;
        end
        else stop:=true;

      key[i+h]:=k;
    end;
  end;
end;

procedure sort_int_vector(key: int_vector);  { SORT a column of numbers}
  var n,part,h,i,j,s: integer;  k,t: integer;  stop: boolean;
begin  { quicksort see knuth Vol 3 }
  n:=high(key);
  part:=1;
  while (twotothe(part+1)-1)<(n div 3) do part:=part+1;
  for s:=part downto 1 do begin
    h:=twotothe(s-1);
    for j:=h+1 to n do begin
      stop:=false;  i:=j-h;  k:=key[j];
      while (i>0) and not stop do
        if k<key[i] then begin
          t:=key[i]; key[i+h]:=t;
          {key[i+h]:=key[i];} i:=i-h;
        end
        else stop:=true;

      key[i+h]:=k;
    end;
  end;
end;

procedure boxlinec(len: integer; line: ansistring);
var i,s: integer; {centered line}
begin
  if len>77 then len:=77;
  s:=(len-length(line)) div 2;
  for i:=1 to (77-len) div 2 do write(' ');
  write('│');
  for i:=1 to s do write(' ');
  write(line);
  for i:=s+length(line)+1 to len do write(' ');
  writeln('│');
end;

procedure boxline(len: integer; line: ansistring);
var i: integer; {not centered line}
begin
  if len>77 then len:=77;
  write('│ ');
  if length(line)>75 then write(line:75) else write(line);
  for i:=length(line)+3 to len do write(' ');
  writeln(' │');
end;

procedure boxtop(len: integer; top: boolean);
const  corner='┌┐└┘';
var i,loc: integer;
begin
  if len>77 then len:=77;
  if top then loc:=1 else loc:=3;
  for i:=1 to (77-len) div 2 do write(' ');
  write(copy(corner,loc,1));  for i:=1 to len do write('─');
  writeln(copy(corner,loc+1,1));
end;

procedure copyright(name,version,years,description: ansistring);
const str1=' Keith W. Kintigh';  str2='All Rights Reserved';
  str3='2014 East Alameda Drive';
  str4='';
  str7='Tempe, Arizona  85282';
var i,maxlen: integer; str5,str6: ansistring;
begin
  {clrscr;}
  for i:=1 to 4 do writeln;
  {gotoxy(1,4);}
  str5:=name+' V'+version;
  str6:='(C) '+years+str1{+', '+str2};
  maxlen:=length(str5);
  if maxlen<length(description) then maxlen:=length(description);
  if maxlen<length(str6) then maxlen:=length(str6);
  maxlen:=maxlen+2;
  boxtop(maxlen,true);
  boxlinec(maxlen,str5);
  boxlinec(maxlen,description);
  boxlinec(maxlen,'');
  boxlinec(maxlen,str6);
  boxlinec(maxlen,str2);
  boxlinec(maxlen,'');
  boxlinec(maxlen,str3);
  boxlinec(maxlen,str7);
  boxtop(maxlen,false);
  for i:=1 to 3 do writeln;
  {gotoxy(1,17);}
end;

procedure displayheading(var fvar: text);
{Display comment lines in first ten lines of a file}
{Whereever used shoud have if filename<>'CON' then displayheading}
var loc,lnecnt: integer;  lne: string[75];  foundcomment: boolean;
begin
  foundcomment:=false;
  lnecnt:=0;
  readln(fvar,lne);
  repeat  {print any of first ten lines with a comment and reset file}
    inc(lnecnt);
    loc:=pos('#',lne);
    if loc>0 then begin
      if not foundcomment then begin
        boxtop(77,true); foundcomment:=true;
      end;
      boxline(77,lne);
    end;
    readln(fvar,lne);
  until eof(fvar) or (lnecnt=10);
  if foundcomment then boxtop(77,false);
  reset(fvar);
end;

{TIME AND DATE FUNCTIONS}

function DateString: ansistring;
begin
  ShortDateFormat:='mm/dd/yy';
  DateString:=DateToStr(Date);
  {getdate(y,m,d,dow);
  begin
    str(y:4,year);
    str(m:2,month);
    str(d:2,day);
  end;
  if month[1]=' ' then month[1]:='0';
  if day[1]=' ' then day[1]:='0';
  year:=copy(year,3,2);
  date := month+'/'+day+'/'+year;}
end;

function ReverseDate(dte: ansistring): ansistring;
var  month,day: string[2];
     year: string[4];
begin
  month:=copy(dte,1,2);
  day:=copy(dte,4,2);
  year:=copy(dte,7,2);
  reversedate:= year+'-'+month+'-'+day;
end;

function timeString: ansistring;
begin
  TimeToStr(Time);
end;

function timesec: real;
begin
  timesec:=frac(time)*24*60*60; {fraction of a 24 hour day}
end;

function readdate(msg: ansistring; dflt: timedatetype): timedatetype;
  {interactively read a date of the form mm/dd/yy - no blanks}
var ok: boolean;
    reply,rslt: timedatetype;
    day,month, year,thisyear: string[2];
    error,value,p1,p2: integer;
begin
  rslt:='';
  thisyear:=copy(DateString,7,2);
  repeat
    ok:=true;
    if dflt<>'' then write(msg,' {',dflt,'} ? ')
    else write(msg,' ? ');
    readln(reply);
    reply:=trim(reply);
    if (reply<>'') or (dflt<>'') then begin
      if length(reply)=0 then reply:=dflt;
      p1:=pos('/',reply);
      if p1<>0 then begin
        reply[p1]:='-';
        p2:=pos('/',reply);
        if p2=0 then begin
          p2:=length(reply)+1;
          reply:=reply+'/'+thisyear;
        end;
        month:=copy(reply,1,p1-1);
        if length(month)=1 then month:='0'+month;
        day:=copy(reply,p1+1,p2-p1-1);
        if length(day)=1 then day:='0'+day;
        year:=copy(reply,p2+1,length(reply)-p2+1);
        if length(year)=1 then year:='0'+year;
        if length(year)=4 then year:=copy(year,3,2);
        val(month,value,error);
        ok:=ok and (error=0) and (value>=1) and (value<=12);
        val(day,value,error);
        ok:=ok and (error=0) and (value>=1) and (value<=31);
        val(year,value,error);
        ok:=ok and (error=0) and (value>0) and (value<99);
        rslt:=month+'/'+day+'/'+year;
      end else ok:=false;
    end;
  until ok;
  readdate:=rslt;
end;

{CSV Procedures}

function gettoken_csv(var ln: ansistring; var p1,error: integer; var tokentype: datatype): ansistring;
{reads a token from the line}
{note in string type embedded blanks are allowwed but may be incompatibe with other programs}
{errors}
    {1 = warning unterminated ", assumed terminated by end of line}
    {2 = " found without preceding separator (blank comma or ^I) separator assumed}
var token: ansistring;
    p2,len,valerror: integer;
    {longintval: longint;}
    extendedval: extended;
    c: ansichar;
    {quote_error: boolean;}
begin
  len:=length(ln);
  token:='';
  error:=0;
  tokentype:=string_type;

  c:=ln[p1];
  case c of
  '"': begin
         if p1=len then begin
           error:=1;  token:='';  p1:=len+1;  {" is last character on line}
         end else begin
           p2:=pos('"',copy(ln,p1+1,len-p1));
           if p2=0 then begin
             error:=1;       {unterminated quote on line}
             token:=copy(ln,p1+1,len-p1); {uses rest of line}
             p1:=len+1;
           end
           else begin {ok, reads to next quote}
             token:=copy(ln,p1+1,p2-2);
             p1:=p1+p2+2;
             while (p1<=len) and (ln[p1]=' ') do inc(p1);
             if (p1<=len) and ((ln[p1]=',') or (ln[p1]=^I)) then inc(p1);
           end;
         end;
       end;
  ',',^I: begin token:=''; inc(p1); end;
  else
    begin
      if p1=len then begin
        token:=c; inc(p1);
      end else begin
        p2:=p1+1;
        while (p2<=len) and (pos(ln[p2],' ,"'^I)=0) do inc(p2);
        token:=copy(ln,p1,p2-p1);
        if (p2<=len)and (ln[p2]='"') then begin
            error:=2; p1:=p2   {" used as separator}
        end else
        if (p2<=len) and (ln[p2]=' ') then begin
          p1:=p2+1;
          while (p1<=len) and (ln[p1]=' ') do inc(p1);
          if (p1<=len) and ((ln[p1]=',') or (ln[p1]=^I)) then inc(p1);
        end else p1:=p2+1;
        token:=trim(token);
      end;
      val(token,extendedval,valerror);  if extendedval=0.0 then;
      if valerror=0 then tokentype:=real_type;
    end;
  end;
  {skip blanks past next separator}
  while (p1<=len) and (ln[p1]=' ') do inc(p1);
  gettoken_csv:=token;
end;

function get_csv(var f: datafile; prompt: ansistring; var filename: ansistring;
quiet: boolean): boolean;
{added quiet parameter, add true for quiet false for verbose in procedure calls}
const maxerrorsbeforehalt=10;
      errormsg: array[1..10] of string=(
        '1:unterminated "',
        '2:unexpected " found w/o preceding separator',
        '3:# of variable labels <> # data values',
        '4:too few data values on line',
        '5:too many data values on line',
        '6: ansistring where real expected',
        '7:unused',
        '8:unused',
        '9:no data lines found',
        '10:internal error'
      );

var token,line: ansistring;
    i,j,p,pos,error,errorcnt,len,obs,nline,lineno,blanklines: integer;
    tokentype: datatype;
    allstring: boolean;
    dskin: text;

  function Nextline: ansistring;
  var s: ansistring;
  begin
    nextline:='';
    repeat   {skip initial blank lines}
      readln(dskin,s);  {get first line}
      inc(lineno);
      s:=trim(s);
    until eof(dskin) or (length(s)>0);
    nextline:=s;
  end;

  procedure csverrorhandler(e,p: integer);
  begin
    if errorcnt<=maxerrorsbeforehalt then begin
      if errorcnt=0 then begin
        writeln;
        writeln(' Line  Pos - #:Error Description');
      end;
      writeln(lineno:5,p:5,' - ',errormsg[e]);
      inc(errorcnt)
    end else begin
      writeln('Program terminated because of too many errors');
      close(dskin);
      closewindow;
      halt;
    end;
  end;

  procedure getcaselabel_csv;
  var i,j,p: integer;
      temp: ansistring;
  begin {ask which string veriable to use and put it in case label}
  with f do begin
    p:=place(tvar);
    if not quiet then begin
      writeln('File Variables:');
      for i:=1 to tvar do begin
        if vindex[i]>0 then write(i:p,'=',padright(vlabel[i],vlabellen))
        else write(i:p,'"',padright(vlabel[i],vlabellen));
        if (i=tvar) or ((i) mod (80 div (vlabellen+p+3))=0) then writeln
        else write('  ');
      end;
      j:=readint('Variable Number of Default Case Label (Reply 0 for case Number)',0,tvar,'0');
    end else j:=0;
    {if j<>0 then} begin
      clabellen:=0;
      if vtype[j]=string_type then
        for i:=1 to nobs do begin
          sdata[i,0]:=truncstring(sdata[i,-vindex[j]],maxclabellen);
          olabel[i]:=sdata[i,0];
          if length(sdata[i,0])>clabellen then clabellen:=length(sdata[i,0]);
        end
      else begin
        p:=readint('Number of Decimal Digits to Use in Label',0,6,'0');
        for i:=1 to nobs do begin
          str(ndata[i,vindex[j]]:nvarplace+P+1:p,temp);
          temp:=trim(temp);
          sdata[i,0]:=truncstring(temp,maxclabellen);
          olabel[i]:=sdata[i,0];
          if length(sdata[i,0])>clabellen then clabellen:=length(sdata[i,0]);
        end;  
      end;
    end;
  end;  end;

begin with f do begin
  {open file}
  filename:=file_prefix(filename)+'.CSV';
  readfile(prompt,dskin,filename);

  {get number of observations}
  nobs:=0; nline:=0; errorcnt:=0;
  while not eof(dskin) do begin
    readln(dskin,line);
    inc(nline);
    line:=trim(line);  {discard empty lines}
    if length(line)>0 then inc(nobs);
  end;
  reset(dskin);

  {retrieve or create variable labels}
  lineno:=0;
  allstring:=true;
  line:=nextline;
  len:=length(line);
  pos:=1;
  tvar:=0;
  vlabellen:=0;
  while (pos<=len) do begin  {decide if these are variable labels}
    token:=gettoken_csv(line,pos,error,tokentype);
    if vlabellen<length(token) then vlabellen:=length(token);
    if tokentype<>string_type then allstring:=false;
    inc(tvar);
  end;

  {allocate tvar label arrays}
  setlength(vlabel,tvar+1);
  setlength(vtype,tvar+1);
  setlength(vindex,tvar+1);

  if allstring then begin {if the first line is labels, extract them}
    if vlabellen>maxvlabellen then vlabellen:=maxvlabellen;
    if vlabellen<4 then vlabellen:=4;
    pos:=1;
    for i:=1 to tvar do begin
      token:=gettoken_csv(line,pos,error,tokentype);
      if error>0 then csverrorhandler(error,pos);
      if length(token)>vlabellen then vlabel[i]:=copy(token,1,vlabellen)
      else vlabel[i]:=token;
    end;
  end else
  begin {create variable labels}
    vlabellen:=place(tvar);
    if vlabellen<4 then vlabellen:=4;
    for i:=1 to tvar do begin
      str(i,token);
      vlabel[i]:='V'+token;
    end;
  end;
  vlabel[0]:='Case';

  blanklines:=nline-nobs;
  if allstring then begin
    dec(nobs); line:=nextline;  {get the next line if we have labels or reinterpret this one}
    if line='' then csverrorhandler(9,0);
  end;

  {infer variable  types from first data line}
  svar:=0; nvar:=0;
  pos:=1; i:=0;  len:=length(line);  vtype[0]:=string_type; vindex[0]:=0;
  while (pos<=len) and (i<tvar) do begin
    token:=gettoken_csv(line,pos,error,tokentype);
    inc(i);
    if tokentype=string_type then begin
      inc(svar);  vtype[i]:=string_type;  vindex[i]:=-svar;
    end else
    begin
      inc(nvar);  vtype[i]:=real_type;  vindex[i]:=nvar;
    end;
  end;

  {check to see if the 1st data line has the same number of values as the label line}
  if not allstring then tvar:=svar+nvar;
  if (tvar<>(svar+nvar)) then begin
    csverrorhandler(3,0);
    for i:=svar+nvar+1 to tvar do begin
      inc(nvar); vtype[i]:=real_type;  vindex[i]:=nvar;
    end;
  end;

  {allocate remainder of arrays}
  setlength(nvlabel,nvar+1);  setlength(nvindex,nvar+1);
  setlength(svlabel,svar+1);  setlength(svindex,svar+1);
  {setlength(clabel,nobs+1);}
  setlength(ndata,nobs+1,nvar+1);  setlength(sdata,nobs+1,svar+1);
  setlength(olabel,nobs+1);

  {finish setting up cross reference indices}
  svindex[0]:=0;  nvindex[0]:=0; vindex[0]:=0; {0 for string not 0 for numeric}
  for i:=1 to tvar do
    if vtype[i]=string_type then begin
      svindex[-vindex[i]]:=i;
      svlabel[-vindex[i]]:=vlabel[i];
    end
    else begin
      nvindex[vindex[i]]:=i;
      nvlabel[vindex[i]]:=vlabel[i];
    end;

  {read, starting with current line, the data file to get data}
  nvarplace:=0;
  obs:=0;
  repeat
    pos:=1; i:=0; len:=length(line);
    inc(obs);
    str(obs,sdata[obs,0]); {create case label;}
    while (pos<=len) and (i<tvar) do begin
      token:=gettoken_csv(line,pos,error,tokentype);
      inc(i);
      if error>0 then csverrorhandler(error,pos);
      if vtype[i]=real_type then begin
        if tokentype=real_type then begin
          val(token,ndata[obs,vindex[i]],error);
          p:=place(ndata[obs,vindex[i]]);
          if p>nvarplace then nvarplace:=p;
        end
        else begin
          csverrorhandler(6,pos);
          ndata[obs,vindex[i]]:=-1;
        end;
      end else sdata[obs,-vindex[i]]:=token;
    end;
    if pos<len then csverrorhandler(5,pos)
    else if i<tvar then begin
      csverrorhandler(4,pos);
      for j:=i+1 to tvar do
        if vtype[j]=string_type then sdata[obs,vindex[j]]:='.'
        else ndata[obs,vindex[j]]:=-1;
    end;
    line:=nextline;
  until line='';
  if obs<>nobs then csverrorhandler(10,0);

  {section to write information about the file}
  writeln('File: ',filename);
  if blanklines>0 then writeln('  ',blanklines,' empty lines ignored');
  if allstring then writeln('  First line interpreted as labels for ',tvar,' variables:');
  writeln(nobs,'  ',' obs & ',tvar,' vars; ',nvar,' numeric & ',svar,' string');
  {end writing section}
  getcaselabel_csv; {get default Case Labels}


  if (filename<>'CON') then close(dskin);
  get_csv:=errorcnt=0;
end; end;

procedure readfile_csv(var df: datafile; prompt: ansistring; var filename: ansistring);
var inputOK,quiet: boolean;
begin
  inputOK:=get_csv(df,prompt,filename,false);
  if not inputOK then begin
    writeln;
    writeln;
    HaltProgram('Input Errors Found, Correct & Rerun');
  end;
end;

procedure readfile_csv_quiet(var df: datafile; prompt: ansistring; var filename: ansistring);
var inputOK: boolean;
begin
  inputOK:=get_csv(df,prompt,filename,true);
  if not inputOK then begin
    writeln;
    writeln;
    HaltProgram('Input Errors Found, Correct & Rerun');
  end;
end;

function ReadFile_csv_Int(var dfi: DataFile_Int; prompt: ansistring; var filename: ansistring): boolean;
var inputOK: boolean;
    i,j: integer;
    dfr: DataFile;
begin
  inputOK:=get_csv(dfr,prompt,filename,false);
  if inputOK then begin
    dfi.tvar:=     dfr.tvar;
    dfi.nvar:=     dfr.nvar;
    dfi.svar:=     dfr.svar;       {number of string variables}
    dfi.nobs:=     dfr.nobs;       {number of observations}
    dfi.clabellen:=dfr.clabellen;  {max length of case labels}
    dfi.vlabellen:=dfr.vlabellen;  {max length of variable labels}
    dfi.nvarplace:=dfr.nvarplace;     {max number of digits to the left of the decimal point, when read}
    dfi.olabel:=   dfr.olabel;
    dfi.vlabel:=   dfr.vlabel;     {all variable labels}
    dfi.svlabel:=  dfr.svlabel;    {string variable labels}
    dfi.nvlabel:=  dfr.nvlabel;    {numeric variable labels}
    setlength(dfi.vtype,dfi.tvar+1);
    for i:=0 to dfi.tvar do dfi.vtype[i]:=dfr.vtype[i];      {type (numeric or string for each varables}
    dfi.sdata:=    dfr.sdata;      {array of string variable values}
    dfi.svindex:=  dfr.svindex;    {svindex[i]=the original variable number of sdata variable i}
    dfi.nvindex:=  dfr.nvindex;    {nvindex[i]=the original variable number of ndata variable i}
    dfi.vindex:=   dfr.vindex;     {vindex[i]:=the index in sdata(0,-) or ndata(+) of original var i}
    setlength(dfi.ndata,dfi.nobs+1,dfi.nvar+1);
    for i:=1 to dfi.nobs do for j:=1 to dfi.nvar do begin
      dfi.ndata[i,j]:=isinteger(dfr.ndata[i,j],inputOK);
      if InputOK then HaltProgram('Non-Integer Data Found in Input File - Numeric Var '+intstrfn(j));
    end;
  end else HaltProgram('Input Errors Found, Correct & Rerun');
end;

procedure datafiletoscreen_csv(var f: datafile);
var i,j,elmtwidth: integer;
begin with f do begin
  writeln;
  if vlabellen<nvarplace+3 then elmtwidth:=nvarplace+3 else elmtwidth:=vlabellen;

  for i:=0 to tvar do begin
    write(vlabel[i]:elmtwidth);
    if (i=tvar) or ((i+1) mod (80 div (elmtwidth+2))=0) then writeln
    else write('  ');
  end;

  for j:=1 to nobs do
    for i:=0 to tvar do begin
      if vtype[i]=real_type then
        write(ndata[j,vindex[i]]:elmtwidth:2)
      else write(truncstring(sdata[j,-vindex[i]],elmtwidth):elmtwidth);
      if (i=tvar) or((i+1) mod (80 div (elmtwidth+2))=0) then writeln
      else write('  ');
    end;
end;  end;

procedure datafiletoprinter_csv(var f:datafile; prompt: ansistring);
var i,j,elmtwidth,ndec,linewid: integer;
    filename: ansistring; dskout: text;
begin with f do begin
  filename:='.TXT';
  writefile(prompt,dskout,filename);
  ndec:=readint('Number of Decimal Places Printed',0,6,'2');
  if vlabellen<(nvarplace+ndec+1) then elmtwidth:=nvarplace+ndec+1 else elmtwidth:=vlabellen;
  linewid:=readint('Output Line Width',elmtwidth,32767,'132');

  for i:=0 to tvar do begin
    write(dskout,vlabel[i]:elmtwidth);
    if (i=tvar) or ((i+1) mod (linewid div (elmtwidth+2))=0) then writeln(dskout)
    else write(dskout,'  ');
  end;

  for j:=1 to nobs do
    for i:=0 to tvar do begin
      if vtype[i]=real_type then
        write(dskout,ndata[j,vindex[i]]:elmtwidth:ndec)
      else write(dskout,truncstring(sdata[j,-vindex[i]],elmtwidth):elmtwidth);
      if (i=tvar) or((i+1) mod (linewid div (elmtwidth+2))=0) then writeln(dskout)
      else write(dskout,'  ');
    end;
  close(dskout);
end;  end;

function csvstring(s: ansistring): ansistring;
var i,len: integer; quote: boolean;
begin
  quote:=false;
  len:=length(s);
  i:=1;
  while not quote and (i<=len) do
  begin
    if s[i]='"' then s[i]:='''';
    quote:=pos(s[i],' ,''')>0;
    inc(i);
  end;
  if quote then csvstring:='"'+s+'"' else csvstring:=s;
end;

function csvreal(v: extended; w,d: integer): ansistring;
var temp: ansistring;
begin
  if frac(v)=0 then str(v:w:0,temp) else str(v:w:d,temp);
  csvreal:=trim(temp);
end;

procedure datafiletofile_csv(var f: datafile; prompt: ansistring);
var i,j,elmtwidth,ndec: integer;
    filename: ansistring; dskout: text;
begin with f do begin
  filename:='.CSV';
  writefile(prompt,dskout,filename);
  ndec:=readint('Number of Decimal Places Retained',0,12,'2');
  if vlabellen<nvarplace+ndec+1 then elmtwidth:=nvarplace+ndec+1 else elmtwidth:=vlabellen;
  for i:=1 to tvar-1 do write(dskout,csvstring(vlabel[i]),',');
  writeln(dskout,csvstring(vlabel[tvar]));

  for j:=1 to nobs do
    for i:=1 to tvar do begin
      if vtype[i]=real_type then
        write(dskout,csvreal(ndata[j,vindex[i]],elmtwidth,ndec))
      else write(dskout,csvstring(sdata[j,-vindex[i]]));
      if (i<tvar) then write(dskout,',') else writeln(dskout);
    end;
  close(dskout);
end;  end;

{MATRIX UTILITIES}

function vector_median(var vect: Real_Vector; elmt:integer): real;
var n,indx: integer;
begin
  n:=high(vect);
  sort_real_vector(vect);
  indx:=n div 2;
  if (n mod 2)=1 then vector_median:=vect[indx+1]
  else vector_median:=(vect[indx]+vect[indx+1])/2.0;
end;

procedure Clear_2D_Int_Array(var table: Int_Array_2D);
var i,j,rows,cols: integer;
begin
  rows:=high(table); cols:=high(table[0]);
  for i:=0 to rows do for j:=0 to cols do table[i,j]:=0;
end;

procedure Fill_Work(var itable: Int_Array_2D; var rtable: Real_Array_2D);
{fill working real table from an integer table}
var i,j,rows,cols: integer;
begin
  rows:=high(itable); cols:=high(itable[0]);
  for i:=0 to rows do for j:=0 to cols do rtable[i,j]:=itable[i,j];
end;

function Int_Vector_Is_0(var vector: Int_Vector): boolean;
var ok: boolean; i,indx: integer;
begin
  ok:=false;
  indx:=high(vector);
  for i:=1 to indx do ok:=ok or (vector[i]=0.0);
  {if ok then writeln('Zero Row or Col');}
  Int_Vector_Is_0:=ok;
end;

procedure Multiply_2D_Real_Array(var table: Real_Array_2D; factor: real);
var i,j,rows,cols: integer;
begin
  rows:=high(table); cols:=high(table[0]);
  for i:=1 to rows do for j:=1 to cols do table[i,j]:=table[i,j]*factor;
end;

procedure Multiply_Real_Row_Or_Col(var table: Real_Array_2D; rowcol: integer;
  rowwise: boolean; factor: real);
var i,j,rows,cols: integer;
begin
  rows:=high(table); cols:=high(table[0]);
  if rowwise then for j:=1 to cols do table[rowcol,j]:=table[rowcol,j]*factor
  else for i:=1 to rows do table[i,rowcol]:=table[i,rowcol]*factor;
end;

procedure Multiply_Real_Vector(var vector: Real_Vector; factor: real);
var j,vals: integer;
begin
  vals:=high(vector);
  for j:=1 to vals do vector[j]:=vector[j]*factor;
end;

function Real_Row_Max(var table: Real_Array_2D; rowno: integer): real;
var i,{rows,}cols: integer;  max: real;
begin
  {rows:=high(table);} cols:=high(table[0]);
  max:=table[rowno,1];
  for i:=2 to cols do
    if table[rowno,i]>max then max:=table[rowno,i];
  Real_Row_Max:=max;
end;

function Real_Column_Max(var table: Real_Array_2D; colno: integer): real;
var i,rows{,cols}: integer;  max: real;
begin
  rows:=high(table); {cols:=high(table[0]);}

  max:=table[1,colno];
  for i:=2 to rows do
    if table[i,colno]>max then max:=table[i,colno];
  Real_Column_Max:=max;
end;

function Real_Vector_Max(var vector: Real_Vector): real;
var i,elmts: integer; max: real;
begin
  elmts:=high(vector);
  max:=vector[1];
  for i:=2 to elmts do if vector[i]>max then max:=vector[i];
  Real_Vector_Max:=max;
end;

function Real_Vector_Min(var vector: Real_Vector): real;
var i,elmts: integer; min: real;
begin
  elmts:=high(vector);
  min:=vector[1];
  for i:=2 to elmts do if vector[i]<min then min:=vector[i];
  Real_Vector_Min:=min;
end;

function Int_Vector_Min(var vector: Int_Vector; elmts: integer): longint;
var i: integer; min: longint;
begin
  min:=vector[1];
  for i:=2 to elmts do if vector[i]<min then min:=vector[i];
  Int_Vector_Min:=min;
end;

function String_Vector_Maxlen(var s: String_Vector): integer;
var i,m,e: integer;
begin
  e:=high(s);
  m:=0;
  for i:=1 to e do if length(s[i])>m then m:=length(s[i]);
  String_Vector_Maxlen:=m;
end;

procedure Margins_2D_Int_Array(var table: Int_Array_2D; var rowsum,colsum: Int_Vector;
  var total: integer);
var i,j,rows,cols: integer;
begin
  rows:=high(table); cols:=high(table[0]);
  total:=0;
  for i:=1 to rows do rowsum[i]:=0;
  for j:=1 to cols do colsum[j]:=0;
  for i:=1 to rows do for j:=1 to cols do begin
    inc(rowsum[i],table[i,j]);
    inc(colsum[j],table[i,j]);
  end;
  for i:=1 to rows do inc(total,rowsum[i]);
end;

procedure Margins_2D_Real_Array(var table: Real_Array_2D; var rowsum,colsum: Real_Vector;
  var total: real);
var i,j,rows,cols: integer;
begin
  rows:=high(table); cols:=high(table[0]);
  total:=0.0;
  for i:=0 to rows do rowsum[i]:=0;
  for j:=0 to cols do colsum[j]:=0;
  for i:=1 to rows do for j:=1 to cols do begin
    rowsum[i]:=rowsum[i]+table[i,j];
    colsum[j]:=colsum[j]+table[i,j];
  end;
  for i:=1 to rows do total:=total+rowsum[i];
end;

procedure pause(var fle: ansistring);
var c: ansichar;
begin
  if fle='CON' then begin
    c:=readchoice('[C]ontinue or [Q]uit','CQ','C');
    if c='Q' then HaltProgram('Program Ends on User Request');
  end;
end;

end.
