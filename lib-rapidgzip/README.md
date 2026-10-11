# lib-rapidgzip

Parallel decompression of ordinary (non-indexed, non-BGZF) gzip files in pure Java, following the
approach of [rapidgzip](https://github.com/mxmlnkn/rapidgzip) (Knespel & Brunst, *Rapidgzip:
Parallel Decompression and Seeking in Gzip Files Using Cache Prefetching*, HPDC 2023).

```java
try (InputStream in = new RapidGzipInputStream(Path.of("reads.fastq.gz"), threads)) {
    // same bytes as new GZIPInputStream(new FileInputStream(...))
}
```

## How it works

1. The compressed file is split into chunks (1 MB by default).
2. For each chunk, a worker thread looks for the first deflate block boundary in it. Every bit
   offset is checked for a valid non-final dynamic Huffman block header: block type, symbol counts,
   a complete precode, valid code length runs and valid literal/length and distance codes with an
   end-of-block code (`BlockFinder`).
3. The worker decodes from there up to the first block boundary in the next chunk without knowing
   the 32 KB of output preceding it. References into that window are output as 16-bit markers;
   once the last 32 KB of output contain no marker, decoding switches to plain bytes
   (`DeflateDecoder`). Gzip member boundaries are handled inside chunks.
4. The reader takes chunks in order. A chunk is only used if it starts exactly where the previous
   chunk ended, which proves its start was a real block boundary. The window for the next chunk is
   resolved from the last 32 KB of this one, and the full marker replacement of the chunk runs on the
   worker pool (`Chunk`, `RapidGzipInputStream`).
5. Whatever speculation can't cover (a false or missed block start) is decoded sequentially with
   the known window.
6. The CRC32 and size in every member trailer are checked, so the output is exactly what a
   sequential decoder produces.

Like `GZIPInputStream`, multiple members are decoded as one stream, the optional header CRC is
checked, and data after the last member which isn't a gzip header is ignored.

It needs random access to the file, so it works on files, not on pipes.

## Build

Java 21+. `./gradlew test` runs the tests, which compare the output with `GZIPInputStream` over
deflate levels and strategies, mixed block types, multiple members, header flags, corrupt and
truncated data, thread counts and chunk sizes.
