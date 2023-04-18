MEMORY
{
    ram (rwx) : ORIGIN = 0x00000000, LENGTH = 8192M /* 128 MB */
}
SECTIONS
{
    . = 0x0;
    .text : {
        *(.text)
    } > ram
    .data : {
        *(.data)
    } > ram
    .sdata : {
	*(.sdata)
    } > ram
    .srodata : {
        *(.srodata)
    } > ram
    .rodata : {
	*(.rodata)
    } > ram
    .riscv.attributes : {
	*(.riscv.attributes)
    } > ram
    .bss : {
        *(.bss)
    } > ram
    .sbss : {
	*(.sbss)
    } > ram
}
