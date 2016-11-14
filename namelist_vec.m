function namelist_vec( fid, name, v )

fprintf( fid, '  %s = ', name );

for i=1:length(v)
    fprintf( fid, '%f', v(i) );

    if ( mod(i, 10) == 0 )
        fprintf( fid, '\n  ' );
    else
        fprintf( fid, '  ' );
    end
end
fprintf( fid, '\n\n' );

end
