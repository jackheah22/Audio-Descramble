module NCO_wrapper(
    input  wire clk,
    output wire [11:0] outsin
);

    // For 7 kHz @ 50 MHz
    localparam PHASE_INC = 32'd601295;  

    simple_dds dds_inst (
        .clk(clk),
        .phase_inc(PHASE_INC),
        .sin_out(outsin)
    );

endmodule