function pno_data = poldata(host_str, port_num)  %poldata('127.0.0.1', 7234);


    % use java Socket and DataInputStream classes
    import java.net.Socket
    import java.io.*    
    
    % connect to the socket and open the data stream
   for attempt = 1:5
        %try
            %fprintf(1, 'Attempt #%d to connect...', attempt);
            socket = Socket(host_str, port_num);
            stream = socket.getInputStream;
            di_stream = DataInputStream(stream);
           % fprintf(1, 'Success\n');
            break;
        % catch
        %     fprintf(1, 'Failure\n');
        %     if ~isempty(socket)
        %         socket.close;
        %     end
        %     if attempt == 5
        %         error('Failed all attempts to connect');
        %     end
        %     pause(1);
        % end % try
    end % for
    
    pno_data = [];
            % get pno data
            pno = zeros(1,6,'single');
            for i = 1:6
                pno(i) = swapbytes(single(di_stream.readFloat) );
            end
            
            % skip crlf
            tmp = zeros(1,2);
            di_stream.read(tmp,0,2);
        
            % create row with data
            pno_row = {pno(1) pno(2) pno(3) pno(4) pno(5) pno(6)};
        
            % append row to pno_data
            pno_data = cat(1,pno_data,pno_row);
      
        socket.close();
end
