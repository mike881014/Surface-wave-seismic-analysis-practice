function [data,dt] = read_svSEG2 (filename)

% read_svSEG2 Read a SEG-2 recorded by StrataVisor(standard SEG-2 format of the Society of Exploration Geophysicist)
%             file from disk.
%
%             [data,dt] = readStrataSEG2 ('filename') reads
%             the file 'filename' and returns the data [m,n] containing n channel of m samples.
%             If no extension is given for the filename, the extension '.dat' is assumed.
%             dt : sampling rate in sec
%
% Written by C-H Lin, Last modified 2007/07/20
% Department of Civil Engineering, National Chiao Tung University


% check argument and filename

if (nargin==0)
   error ('SEG2LOAD requires at least a filename as an argument !');
end;

if (isstr (filename)~=1)
   error ('Argument is not a filename !');
end;

if (isempty (findstr (filename,'.'))==1)
   filename =[filename,'.dat'];
end;

fid=fopen (filename,'rb','ieee-le');
if (fid ==-1)
    error (['Error opening ',filename,' for input !']);
end;

% check for SEG-2 file type
% first 2 bytes equal '3a55h' (14933) for PC/Windows
% first 2 bytes equal '553ah' (21818) for UNIX

  fileType=fread (fid,1,'short');
  if (fileType ~= 14933)
    fclose (fid);
    error ('Not a SEG-2 file !');
  end;

% read file descriptor block

  revNumber          = fread (fid,1,'short');
  sizeOfTracePointer = fread (fid,1,'ushort');
  nbOfTraces         = fread (fid,1,'ushort');
  sizeOfST           = fread (fid,1,'uchar');
  firstST            = fread (fid,1,'char');
  secondST           = fread (fid,1,'char');
  sizeOfLT           = fread (fid,1,'uchar');
  firstLT            = fread (fid,1,'char');
  secondLT           = fread (fid,1,'char');
  reserved           = fread (fid,18,'uchar');
  tracePointers      = fread (fid,nbOfTraces,'ulong');  

% read free strings 
  
  fseek (fid,32+sizeOfTracePointer,'bof');
  offset = fread (fid,1,'ushort');
  while (offset > 0)
    freeString = char (fread (fid,offset-2,'char'))';
    offset = fread (fid,1,'ushort');
  end;

% read traces

  for i=1:nbOfTraces,
    fseek (fid,tracePointers (i),'bof');
    traceId     = fread (fid,1,'ushort');
    sizeOfBlock = fread (fid,1,'ushort');
    sizeOfData  = fread (fid,1,'ulong');
    nbOfSamples = fread (fid,1,'ulong');
    dataCode    = fread (fid,1,'uchar');
    reserved    = fread (fid,19,'uchar');
    freeString  = char(fread (fid,sizeOfBlock-32,'char'))';
    data (:,i)  = fread (fid,nbOfSamples,'float32');
  end;

  %freeString
  %freeString((findstr(freeString,'SAMPLE_INTERVAL')+length('SAMPLE_INTERVAL')):(findstr(freeString,'SHOT_SEQUENCE_NUMBER')-1))
%if (findstr (freeString,'SAMPLE_INTERVAL') > 0)
    dt = str2num (freeString((findstr(freeString,'SAMPLE_INTERVAL')+length('SAMPLE_INTERVAL')):(findstr(freeString,'SHOT_SEQUENCE_NUMBER')-4)))
    %end;

fclose (fid);
