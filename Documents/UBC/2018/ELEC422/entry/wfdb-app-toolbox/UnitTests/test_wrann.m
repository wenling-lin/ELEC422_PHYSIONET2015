function [tests,pass,perf]=test_wrann(varargin)

inputs={'verbose'};
verbose=0;
for n=1:nargin
    if(~isempty(varargin{n}))
        eval([inputs{n} '=varargin{n};'])
    end
end


%Test the examples 
test1_str=[ '[ann,type,subtype,chan,num]=rdann(''challenge/2013/set-a/a01'',''fqrs'');' ...
      'wrann(''challenge/2013/set-a/a01'',''test'',ann,type,subtype,chan,num);' ...
      '[ann,type,subtype,chan,num]=rdann(''challenge/2013/set-a/a01'',''fqrs'');' ...
      'wrann(''challenge/2013/set-a/a01'',''test'',ann,type,subtype,chan,num);' ...
      '[ann2,type2,subtype2,chan2,num2]=rdann(''challenge/2013/set-a/a01'',''test'',[],[],1);' ...
      'err=sum(ann ~= ann2);'];

test2_str=[ '[ann,type,subtype,chan,num]=rdann(''mitdb/100'',''atr'');' ...
	  'wrann(''mitdb/100'',''test'',ann,type,subtype,chan,num);'];

%Test that comments are being written
test3_str=['[ann,type,subtype,chan,num,comments]=rdann(''afdb/04015'',''atr'');' ...
            'wrann(''afdb/04015'',''test'',ann,type,subtype,chan,num,comments);' ...
            '[ann2,type2,subtype2,chan2,num2, comments2]=rdann(''afdb/04015'',''test'');'...
            'equal=strcmp([comments{:}],[comments2{:}]);' ...
            'if(equal==0);error(''comments are incorrect'');end;' ];

test_string={test1_str,test2_str, test3_str};
clean_up={['delete([pwd filesep ''challenge'' filesep ''2013'' filesep ''set-a'' filesep ''*'']);' ...
          'rmdir([pwd filesep ''challenge''],''s'');'], ...
          ['delete([pwd filesep ''mitdb'' filesep ''100'' filesep ''*'']);' ...
          'rmdir([pwd filesep ''mitdb''],''s'');'] ...
          ['delete([pwd filesep ''afdb'' filesep ''*'']);' ...
          'rmdir([pwd filesep ''afdb''],''s'');'], ...
          };
[tests,pass,perf]=test_wrapper(test_string,clean_up,verbose);

 
 