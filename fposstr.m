function found = posfile( fid, searchstr )


found = false;
while ~feof( fid )

    tline = fgetl( fid );

    if ischar( tline ) && ~isempty( strfind( tline, searchstr ))
        found = true;
        break;
    end

end

return
end
