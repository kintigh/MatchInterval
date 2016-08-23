program MatchInterval;

{$APPTYPE CONSOLE}  {Delphi & Lazarus}
{$DEFINE LAZ}

uses
  {Delphi & Lazarus:}  SysUtils, strutils,
  {IFDEF LAZ KWKStdXEL {ELSE}
  {Delphi} KWKStdXE {ENDIF} ;

const ver='2.2';  {2.2 eliminat echeck on analysis interval start/end overalapping fixed interval}
                  {2.1 MatchXErad3 debugged from MatchXErad2 (called 1.1); Renamed to MatchInterval2_1}
                  {2.0 Matchinterval seems to be earlier than MatchXErad2, so ignored}
                  {1.1 MatchXErad2 directe ancestor of 2.1, adds auto comput date selection}
                  {1.0 2013-07-08 Start}
      years='2013-2016' ;
      debug=false;
      max_rule=3;
      gt1rule=false;
      minmoviningintervalrule=false;


var  df: datafile;
     Last_Fixed_Interval,Last_Moving_Interval,Gap_Intervals,First_Fixed_Interval,First_Moving_Interval,
       rule,earliestF,latestF,earliestM,latestM,computation_start,computation_end,analysis_start,analysis_end,
       Random_Runs,minmovingintervalstart,maxmovingintervalstart,minmovingintervallength,minoverlap,
       totalextremeyears,impacttimelag,negfixedstartuncertainty,posfixedstartuncertainty:integer;
     Moving_Interval_Start,Moving_Interval_End,Moving_Interval_Length,
       Save_Moving_Interval_Start,Save_Moving_Interval_End,Save_Moving_Interval_Length,
       Fixed_Interval_Start,Fixed_Interval_End,Gap_Length,Save_Gap_Length: array of integer;
     Moving_Interval_Random: array of real;
     No_Gap_At_Beginning,print_fixed,Allow_GT_1_Match,AdjustComputationDates: boolean;
     listfile,filein1,filein2: Ansistring;
     list: text;

Procedure allocate_arrays(size:integer);
begin
  setlength(Moving_Interval_Start,size);
  setlength(Moving_Interval_End, size);
  setlength(Save_Moving_Interval_Start,size);
  setlength(Save_Moving_Interval_End,  size);
  setlength(Moving_Interval_Length,size);
  setlength(Save_Moving_Interval_Length,size);
  setlength(Gap_Length,size+2);
  setlength(Save_Gap_Length,size+2);
  setlength(Moving_Interval_Random,size+2);
end;

Procedure process_randomized_events;
var i,nvar: integer;
begin
  Last_Moving_Interval:=df.nobs;
  nvar:=df.nvar;

  if (nvar <> 2) then
    HaltProgram('Randomized Event: Invalid Number of Input Variables.');

  allocate_arrays(Last_Moving_Interval+1);

  {setlength(Moving_Interval_Start,Last_Moving_Interval+1);
  setlength(Moving_Interval_End,Last_Moving_Interval+1);
  setlength(Moving_Interval_Length,Last_Moving_Interval+1);
  setlength(Save_Moving_Interval_Length,Last_Moving_Interval+1);
  setlength(Gap_Length,Last_Moving_Interval+3);
  setlength(Save_Gap_Length,Last_Moving_Interval+3);
  setlength(Moving_Interval_Random,Last_Moving_Interval+3); }

  earliestM:=maxint;
  latestM:=-maxint;

  for i:=1 to Last_Moving_Interval do begin
    Moving_Interval_Start[i]:=round(df.ndata[i,1]);
    Moving_Interval_End[i]:=round(df.ndata[i,2]);

    if Moving_Interval_Start[i]<earliestM then earliestM:=Moving_Interval_Start[i];
    if Moving_Interval_End[i]>latestM then latestM:=Moving_Interval_End[i];

    if debug then writeln(list,'Moving',i:4,Moving_Interval_Start[i]:6,Moving_Interval_End[i]:6);

    if Moving_Interval_Start[i]>Moving_Interval_End[i] then
      HaltProgram('Randomized Event Start Date After End Date');
    {program does not check for overlapping intervals}
  end;
end;

Procedure process_fixed_events;
var i,nvar: integer;
begin
  Last_Fixed_Interval:=df.nobs;
  nvar:=df.nvar;

  if (nvar <> 2) then
    HaltProgram('Fixed Event: Invalid Number of Input Variables.');

  setlength(Fixed_Interval_Start,Last_Fixed_Interval+1);
  setlength(Fixed_Interval_End,Last_Fixed_Interval+1);
  earliestF:=maxint;
  latestF:=-maxint;

  for i:=1 to Last_Fixed_Interval do begin
    Fixed_Interval_Start[i]:=round(df.ndata[i,1]);
    Fixed_Interval_End[i]:=round(df.ndata[i,2]);

    if Fixed_Interval_Start[i]<earliestF then earliestF:=Fixed_Interval_Start[i];
    if Fixed_Interval_End[i]>latestF then latestF:=Fixed_Interval_End[i];

    if debug then writeln(list,'Fixed',i:4,Fixed_Interval_Start[i]:6,Fixed_Interval_End[i]:6);
    if Fixed_Interval_Start[i]>Fixed_Interval_End[i] then
      HaltProgram('Randomized Event Start Date After End Date');

    {program does not check for overlapping intervals}
  end;
end;

procedure prestart;
begin
   copyright('MatchIntervals',ver,years,
    'Monte Carlo Assessment of Correspondence Between 2 Sets of Intervals');
  if debug then listfile:='CON' else listfile:='.TXT';
  writefile('Results File',list,listfile);
  writeln(list,'Program MatchIntervals v',ver,' - Keith Kintigh - 20',reversedate(datestring),' ',TimeToStr(time))  ;
  writeln(list);
end;

procedure get_data;
Begin
  writeln;
  {Get fixed Inverval Data - eg social transformations}
  if debug then filein1:='TestF' else filein1:='';
  readfile_csv_quiet(df,'Fixed (Transformations) Interval File (Start, End)',filein1);
  process_fixed_events;
  writeln;
  if debug then filein2:='TestM' else filein2:='';      {fix}
  readfile_csv_quiet(df,'Moving (Climate) Interval File (Start, End)',filein2);
  process_randomized_events;
  writeln;
  writeln('Start of Earliest Fixed Interval ',earliestF:6);
  writeln('End of Latest Fixed Interval     ',latestF:6);
  writeln('Start of Earliest Moving Interval',earliestM:6);
  writeln('End of Latest Moving Interval    ',latestM:6);
End;

procedure calculate_gaps;
var i: integer;
begin
  No_Gap_At_Beginning:=(Moving_Interval_Start[First_Moving_Interval]=computation_start);  {Start sequence with target event}

  Gap_Intervals:=1;
  if not No_Gap_At_Beginning then begin
    Gap_Length[1]:=Moving_Interval_Start[1]-computation_start;  inc(Gap_Intervals);
  end;

  for i:=1 to Last_Moving_Interval do begin
    Moving_Interval_Length[i]:=Moving_Interval_End[i]-Moving_Interval_Start[i]+1;
    if i<Last_Moving_Interval then begin
      Gap_Length[Gap_Intervals]:=Moving_Interval_Start[i+1]-Moving_Interval_End[i]-1;
      inc(Gap_Intervals);
    end;
  end;
  if Moving_Interval_End[Last_Moving_Interval]<>computation_end then
    Gap_Length[Gap_Intervals]:=computation_end-Moving_Interval_End[Last_Moving_Interval]
  else dec(Gap_Intervals);
end;

procedure initialize;
var i,original_last: integer;  s: string;   error: boolean;
begin
  error:=true;
  while error do begin
    writeln;
    str(earliestM,s);
    analysis_start:=readint('Analysis Interval Start Date    ',-40000,2100,s);
    str(latestM,s);
    analysis_end:=readint(  'Analysis Interval End Date      ', analysis_start+1,2101,s);
    error:=false;
    {for i:=1 to last_fixed_interval do begin
       error:=error or ((analysis_start>=fixed_interval_start[i]) and (analysis_start<=fixed_interval_end[i])) or
          ((analysis_end>=fixed_interval_start[i]) and (analysis_end<=fixed_interval_end[i]));
       if error then writeln('Warning: the analysis interval start or end is during a fixed interval (T',i,')');
    end;}
  end;
  computation_start:=analysis_start;
  computation_end:=analysis_end;

  if (analysis_start<earliestF) or (analysis_end>latestF) then AdjustComputationDates:=true else
    AdjustComputationDates:=Readbool('Adjust Analysis Intervals to Start and End of Fixed Intervals','T');

  {ignore moving intervals before and after start and end}
  First_Moving_Interval:=1;
  while Moving_Interval_End[First_Moving_Interval]<computation_start do inc(First_Moving_Interval);      {ignore interval}
  if Moving_Interval_Start[First_Moving_Interval]<computation_start then begin {analysis interval starts during a moving interval}
    if AdjustComputationDates then computation_start:=Moving_Interval_Start[First_Moving_Interval]   {adjust global to beginning of moving interval}
    else Moving_Interval_Start[First_Moving_Interval]:=computation_start;                  {truncate moving interval}
  end else
  if computation_start<Moving_Interval_Start[First_Moving_Interval] then begin {analysis interval starts during a gap}
    if AdjustComputationDates then begin
      if (First_Moving_Interval>1) then computation_start:=Moving_Interval_End[First_Moving_Interval-1]+1 {start at beginning of previous gap}
      else computation_start:=Moving_Interval_Start[First_Moving_Interval];                               {Start at beginning of next moving interval}
    end
  end;
  original_Last:=Last_Moving_Interval;
  while Moving_Interval_Start[Last_Moving_Interval]>computation_end do dec(Last_Moving_Interval);        {ignore interval}
  if Moving_Interval_End[Last_Moving_Interval]>computation_end then begin
    if AdjustComputationDates then computation_end:=Moving_Interval_End[Last_Moving_Interval] {analysis ends during a moving interval}
    else Moving_Interval_End[Last_Moving_Interval]:=computation_end;
  end else
  if computation_end>Moving_Interval_End[Last_Moving_Interval] then begin {analysis ends during a gap}
    if AdjustComputationDates then begin
      if (Last_Moving_Interval<Original_Last) then computation_end:=Moving_Interval_Start[Last_Moving_Interval+1]-1 {start at end of the next ga}
      else computation_end:=Moving_Interval_End[Last_Moving_Interval];                               {Start at beginning of next moving interval}
    end
  end;

  {Discard outside intervals}
  if First_Moving_Interval<>1 then begin
    for i:=First_Moving_Interval to Last_Moving_Interval do begin
      Moving_Interval_Start[i-First_Moving_Interval+1]:=Moving_Interval_Start[i];
      Moving_Interval_End[i-First_Moving_Interval+1]:=Moving_Interval_End[i];
    end;
    Last_Moving_Interval:=Last_Moving_Interval-First_Moving_Interval+1;
    First_Moving_Interval:=1;
  end;

  {ignore fixed intervals before and after start and end}
  First_Fixed_Interval:=1;
  while Fixed_Interval_End[First_Fixed_Interval]<computation_start do inc(First_Fixed_Interval);      {ignore interval}
  if Fixed_Interval_Start[First_Fixed_Interval]<computation_start then
   Fixed_Interval_Start[First_Fixed_Interval]:=computation_start;  {truncate}
  while Fixed_Interval_Start[Last_Fixed_Interval]>computation_end do dec(Last_Fixed_Interval);              {ignnoe interval}
  if Fixed_Interval_End[Last_Fixed_Interval]>computation_end then
    Fixed_Interval_End[Last_Moving_Interval]:=computation_end;     {truncate}

  {Discard outside intervals}
  if First_Fixed_Interval<>1 then begin
    for i:=First_Fixed_Interval to Last_Fixed_Interval do begin
      Fixed_Interval_Start[i-First_Fixed_Interval+1]:=Fixed_Interval_Start[i];
      Fixed_Interval_End[i-First_Fixed_Interval+1]:=Fixed_Interval_End[i];
    end;
    Last_Fixed_Interval:=Last_Fixed_Interval-First_Fixed_Interval+1;
    First_Fixed_Interval:=1;
  end;

  totalextremeyears:=0;
  for i:=1 to last_moving_interval do
    inc(totalextremeyears,moving_interval_end[i]-moving_interval_start[i]+1);

  calculate_gaps;

  if AdjustComputationDates then writeln('(Adjusted) Analysis Interval: ',computation_start,' to ',computation_end);

  writeln;
  if minmoviningintervalrule then minmovingintervallength:=readint('Minimum Moving Interval Length to Consider',0,maxint,'1')
  else minmovingintervallength:=1;

  if gt1rule then Allow_GT_1_Match:=readbool('Allow >1 Moving Intervals to Match a Given Fixed Interval','F')
  else Allow_GT_1_Match:=false;

  writeln;
  Random_Runs:=readint('Number of Random Runs',1,10000000,'1000000');
  setrandom;
  writeln;
end;
 

Procedure get_rule;
begin
  writeln;
  rule:=pos(readchoice('Match: [O]verlap, Moving [L]ag Interval; Moving Lag/Fixed [U]ncertainty','OLU','U'),'OLU');

  if rule=1 then minoverlap:=readint('  Minimum Overlap Required for Match',0,maxint,'1');
  if rule=2 then begin
    writeln('  Fixed Interval must start at least [Min] and no more than [Max] years');
    writeln('  after the start of the Moving Interval');

    minmovingintervalstart:=readint('  Min Years after Start of Moving Interval for Fixed Interval Start ',
      -maxint,maxint,'5');
    maxmovingintervalstart:=readint('  Max Years after Start of Moving Interval for Fixed Interval Start',
      minmovingintervalstart,maxint,'5');
  end;
  if rule=3 then begin

    ImpactTimeLag:=readint('Impact Time Lag from Beginning of Moving Interval',
      -maxint,maxint,'5');
    if readbool ('  Symmetric Uncertainty around Fixed Interval Start','T') then begin
      posfixedstartuncertainty:=readint('  Fixed Start Uncertainty +/-',0,maxint,'5');
      negfixedstartuncertainty:=-posfixedstartuncertainty;
    end
    else begin
      negfixedstartuncertainty:=readint('  - Fixed Start Uncertainty: Years Before Dated Start',0,maxint,'5');
      posfixedstartuncertainty:=-readint('  + Fixed Start Uncertainty: Years After Dated Start',0,maxint,'5');
    end;
  end;
end;

procedure data_list(full: boolean);
var i,g: integer;
begin
  writeln(list,'================================');
  if full then begin
    writeln(list,'All Interval Data');
    writeln(list);
    writeln(list,'Fixed Interval File:  ',filein1);
    writeln(list,'Moving Interval File: ',filein2);
  end
  else Begin
    writeln(list,'Intervals Evaluated');
    writeln(list,'  Note: Intervals before Analysis Start and End Dates are Ignored');
  End;
  print_fixed:=true;

  g:=1;
  if print_fixed then begin
    writeln(list);
    writeln(list,'      Fixed Intervals');
    writeln(list,'   Start     End   Years');

    for i:=1 to Last_Fixed_Interval do Begin
      writeln(list,Fixed_Interval_Start[i]:8,Fixed_Interval_End[i]:8,
        Fixed_Interval_End[i]-Fixed_Interval_Start[i]+1:8);
    end;

    print_fixed:=false;
  end;
  writeln(list);
  write(list,'      Moving Intervals');
  if not full then writeln(list,'      Next') else writeln(list);
  write(list,'   Start     End   Years');
  if not full then writeln(list,'     Gap') else writeln(list);
  if (not full) and (not no_gap_at_beginning) then begin
    writeln(list,' ':24,Gap_Length[g]:8);
    inc(g);
  end;
  for i:=1 to Last_Moving_Interval do begin
    write(list,Moving_Interval_Start[i]:8,Moving_Interval_End[i]:8,Moving_Interval_End[i]-Moving_Interval_Start[i]+1:8);
    if (not full) and (g<=gap_intervals) then writeln(list,Gap_Length[g]:8)else writeln(list);
    inc(g);
  end;
  writeln(list,'================================');
  {writeln(list,'Gap Interval Length');
  for i:=1 to Gap_Intervals do writeln(list,Gap_Length[i]);}
end;

Procedure print_file_parms;
begin
  writeln(list,'Start of Earliest Fixed Interval ',earliestF:6);
  writeln(list,'End of Latest Fixed Interval     ',latestF:6);
  writeln(list,'Start of Earliest Moving Interval',earliestM:6);
  writeln(list,'End of Latest Moving Interval    ',latestM:6);
  writeln(list);
  writeln(list,'Analysis Start: ',analysis_start,'   End: ',analysis_end,
    '   Years: ',analysis_end-analysis_start+1);
  writeln(list);

  writeln(list,'Computation Start: ',computation_start,'   Computation End: ',computation_end,
    '   Years: ',computation_end-computation_start+1);
  if AdjustComputationDates then begin
    writeln(list,'  Computation Start is the start of the moving interval or gap including ',analysis_start);
    Writeln(list,'  Computation End is the end of the moving interval or gap including ',analysis_end);
  end;
  writeln(list);
  writeln(list,'Years in Moving Intervals: ',totalextremeyears:4,' or ',
    100.0*totalextremeyears/(computation_end-computation_start+1):5:1,'% of Computation Interval');
  if No_Gap_At_Beginning then writeln(list,'Don''t start with Gap')
  else writeln(list,'Start with Gap');

  writeln(list,'Moving intervals must be at least ',minmovingintervallength,' years long to be considered');
  if allow_gt_1_match then writeln(list,'Count >1 Match/Fixed Interval')
  else writeln(list,'Only Count 1 Match/Fixed Interval');

  writeln(list,'Random Runs: ',Random_Runs:8);
  writeln(list,'Random Number Generator Seed: ',randseed);

  writeln(list);
end;

procedure print_parms;
begin
  writeln(list);

  {writeln(list,'Matching Rule ',rule:6);}
  case rule of
    1:  begin
      writeln(list,'Rule 1 - Overlap: ',minoverlap);
      writeln(list,'     Overlapping Intervals');
    end;
    2: begin
      writeln(list,'Rule 2 - Moving Starts ',minmovingintervalstart,' - ',minmovingintervalstart, ' after Fixed Start');
      writeln(list,'     Moving Interval Start Date Relative to Fixed Interval Start');
      writeln(list,'     Fixed Interval must start at least ',minmovingintervalstart,' years and');
      writeln(list,'     and no more than ',minmovingintervalstart,' years after the start of the Moving Interval');
    end;
    3: begin
      writeln(list,'Rule 3 - Lag: ',impacttimelag,'  Uncertainty ', negfixedstartuncertainty,' : ',posfixedstartuncertainty);
      writeln(list,'     Impact Time Lag & Fixed Interval Uncertainty');
      writeln(list,'     Moving Interval Impact Time Lag: ',impacttimelag,' from Start of Moving Interval');
      writeln(list,'     Impact is Up To ',-negfixedstartuncertainty,' Years Before Dated Fixed Interval Start Date and');
      writeln(list,'     Impact is Before ',posfixedstartuncertainty,' Years After Dated Fixed Interval Start Date');
    end;
  end;
  writeln(list);
end;

procedure restore_actual;
var i: integer;
begin
  for i:=1 to Last_Moving_Interval do begin
    Moving_Interval_Start[i] := Save_Moving_Interval_Start[i];
    Moving_Interval_End[i]   := Save_Moving_Interval_End[i];
    Moving_Interval_Length[i]:= Save_Moving_Interval_Length[i];
  end;
  for i:=1 to Gap_Intervals do Gap_Length[i]:=Save_Gap_Length[i];
end;

procedure save_actual;
var i: integer;
begin
  for i:=1 to Last_Moving_Interval do begin
    Save_Moving_Interval_Start[i] := Moving_Interval_Start[i];
    Save_Moving_Interval_End[i]   := Moving_Interval_End[i];
    Save_Moving_Interval_Length[i]:= Moving_Interval_Length[i];
  end;
  for i:=1 to Gap_Intervals do Save_Gap_Length[i]:=Gap_Length[i];
end;

procedure sort_list(var key: array of real; var rec: array of integer; n,part: integer);  { SORT a column of numbers}
var h,i,j,r,s: integer;  k: real;  stop: boolean;
begin  { quicksort see knuth Vol 3 }
  for s:=part downto 1 do begin
    h:=twotothe(s-1);
    for j:=h+1 to n do begin
      stop:=false;  i:=j-h;  k:=key[j];  r:=rec[j];
      while (i>0) and not stop do
        if k<key[i] then begin
          key[i+h]:=key[i];  rec[i+h]:=rec[i];  i:=i-h;
        end
        else stop:=true;
      key[i+h]:=k;  rec[i+h]:=r;
    end;
  end;
end;

function overlap(startdate1,enddate1,startdate2,enddate2: integer):integer;
{compute the integer overlap between two intervals}
begin
  overlap:= i_maximum(0, i_minimum(EndDate1, EndDate2) - i_maximum(StartDate1, StartDate2) +1 );
end;

function within(ints,inte,target: integer): boolean;
{see if target date is within interval inclusive of end points}
begin
  within:=(ints<=target) and (target<=inte);
end;

procedure displaymatch(fs,fe,ms,me:integer; dupe:boolean);
begin
  write(list,'Match: Fixed  ',fs:4,'-',fe:4,'  Moving ',ms:4,'-',me:4);
  if dupe then writeln(list,' (Duplicate Match)') else writeln(list);
end;

function match(run: integer): integer;
var i,j,nmatch:integer; fixedmatch,thesematch: boolean;
begin
  nmatch:=0;
  thesematch:=false;
  for j:=1 to Last_Fixed_Interval do begin
    fixedmatch:=false;
    for i:=1 to Last_Moving_Interval do begin
      if (Moving_Interval_End[i]-Moving_Interval_Start[i]+1)>=minmovingintervallength then
      case rule of
        1: thesematch:=overlap(Fixed_Interval_Start[j],Fixed_Interval_End[j],
             Moving_Interval_Start[i],Moving_Interval_End[i]) >= minoverlap;
        2: thesematch:=within(Moving_Interval_Start[i]+minmovingintervalstart,
             Moving_Interval_Start[i]+maxmovingintervalstart,Fixed_Interval_Start[j]);
        3: thesematch:=within(Fixed_Interval_Start[j]+negfixedstartuncertainty,
             Fixed_Interval_Start[j]+posfixedstartuncertainty,Moving_Interval_Start[i]+ImpactTimeLag);
      end;
      if thesematch then begin
        if (not fixedmatch) or allow_gt_1_match then inc(nmatch);
        if debug or (run=0) then displaymatch(Fixed_Interval_Start[j],Fixed_Interval_End[j],
          Moving_Interval_Start[i],Moving_Interval_End[i],fixedmatch);
        fixedmatch:=true;
      end;
    end;
  end;
  match:=nmatch;
end;


procedure analysis;
var i,r,Trial_Match,part,Gap_Index,offset,actual_matches,total_matches,
    matches_ge_actual: integer;
    prob: real;
begin
  actual_matches:=match(0);    {sees how many matches in actual data}
  matches_ge_actual:=0;     {counter for random runs}
  total_matches:=0;
  writeln(list);

  part:=1;   {this is setup for the quicksort}
  while (twotothe(part+1)-1)<(Gap_Intervals div 3) do part:=part+1;

  for r:=1 to Random_Runs do begin   {random runs here}
    if debug then begin
      writeln(list);
      writeln(list,'---------------------------');
      writeln(list,'Run ',r);
    end;
    for i:=1 to Last_Moving_Interval do begin
      Moving_Interval_Random[i]:=random;    {randomize extreme events}
    end;
    sort_list(Moving_Interval_Random,Moving_Interval_Length,Last_Moving_Interval,part);
    for i:=1 to Gap_Intervals do Moving_Interval_Random[i]:=random;     {randomize normal events}
    sort_list(Moving_Interval_Random,Gap_Length,Gap_Intervals,part);

    if debug then begin
      Writeln(List,'Sort');
      for i:=1 to Last_Moving_Interval do writeln(List,'  Moving ',Moving_Interval_Length[i]:3);
      for i:=1 to Gap_intervals do    writeln(list,'  Gaps   ',Gap_Length[i]:3);
    end;

    if No_Gap_At_Beginning then begin offset:=0; Gap_Index:=1; end  {construct randomized sequence}
    else begin offset:=Gap_Length[1]; Gap_Index:=2; end;    {of extreme and normal intervals}
    for i:=1 to Last_Moving_Interval do begin
        Moving_Interval_Start[i]:=computation_start+offset;           {we only keep track of extremes}
        Moving_Interval_End[i]:=Moving_Interval_Start[i]+Moving_Interval_Length[i]-1;
        offset:=offset+Moving_Interval_Length[i]+Gap_Length[Gap_Index];
        inc(Gap_Index);
    end;

    if debug then data_list(false);

    Trial_Match:=match(r);    {See if there is a match}

    {debugging code}
    {if trial_match>1 then begin
      writeln(list,'Trial_Match>1');
      data_list;
    end;}

    if debug then writeln(list,'Matches ',Trial_Match);

    total_matches:=total_matches+Trial_Match ;
    if Trial_Match>=actual_matches then inc(matches_ge_actual);

  end;

  prob:=matches_ge_actual/Random_Runs;

  writeln(list,'Observed    Random Probability');
  writeln(list,' Matches   Matches Rand.>=Obs.');
  writeln(list,actual_matches:8,total_matches/Random_Runs:10:4,prob:12:4);
  writeln(list);

  writeln;
  writeln('Observed    Random Probability');
  writeln(' Matches   Average Random>=Obs');
  writeln(actual_matches:8,total_matches/Random_Runs:10:4,prob:12:4);
  writeln;
  writeln(list,'--------------------------------');
end;

procedure wrapup;
begin
  writeln(list,'Program End');
  close(list);
  writeln('Program End');
end;


begin
  {main program block}
  prestart;
  repeat
    get_data;
    data_list(true);
    initialize;
    print_file_parms;
    data_list(false);
    save_actual;
    repeat
      get_rule;
      print_parms;
      restore_actual;
      analysis;
    until not readbool('Run Again with Different Decision Rules','F');
  until not readbool('Run Again with Different Data','F');
  wrapup;
  CloseWindow;
end.

