function logger(to_log)
% Log string 'string'
	global g_logfile;
	string = sprintf('[%s] %s\n',datestr(now,0),to_log);
	fprintf(string);
   	fprintf(g_logfile,string);
end
