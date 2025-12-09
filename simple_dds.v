module simple_dds (
    input  wire        clk,
    input  wire [31:0] phase_inc,
    output reg  [11:0] sin_out
);

    reg [31:0] phase = 0;

    always @(posedge clk) begin
        phase <= phase + phase_inc;
        sin_out <= sine_lut[phase[31:20]];  
        // phase[31:20] = top 12 bits = 4096-entry table index
    end

    // 4096-entry sine lookup table (12-bit samples)
    // You generate this once with Python/Matlab.
    reg [11:0] sine_lut [0:4095];
    initial $readmemh("sine4096.hex", sine_lut);

endmodule