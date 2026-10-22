# 0004. Retargetable Computation Stack

- **Status:** Accepted; amended by [0008](0008-independent-libraries.md): the stack now lives in calx
- **Date:** 2026-10-21
- **Decided by:** the user set the direction: a conservative, TVM-like separation, CPU first, open to future devices. Claude designed the components.

## 1. Context

The user mainly computes on CPUs, but wants computations to be deployable to many devices, including devices that do not yet exist. Hardware changes faster than scientific code. If a design names today's devices in its core abstractions, it will need rewriting with every new device.

## 2. Decision

Build `core/compute/` as a separated stack in the manner of TVM, Halide and MLIR ([architecture §4](../../docs/architecture.md)). It has six components:

1. **Representation** (`ir/`): device-neutral, versioned and serialisable.
2. **Schedule** (`schedule/`): how one computation is mapped onto hardware, kept separate from what it computes.
3. **Compiler** (`compiler/`): passes and lowering.
4. **Runtime** (`runtime/`): a minimal device API with a C ABI and DLPack arrays.
5. **Scheduler** (`scheduler/`): task graphs across devices.
6. **Backends** (`backends/`): plug-ins, each with a target description, a code generator and a device-API implementation.

Further rules:

- **The CPU is the first backend.** Its reference interpreter is the correctness oracle for all other backends.
- **The project owns these interfaces.** External systems, such as LLVM, MLIR, TVM, Kokkos and SYCL implementations, may be used inside backends. They never become the contract.
- **Scalar algorithms stay ordinary C++.** Only kernels that need several devices go through the IR.
- **Precision is declared per computation**, and is never silently reduced.

## 3. Alternatives Considered

- **Kokkos or SYCL as the abstraction.** Rejected as the contract, because it ties the representation to C++ templates and to the devices supported today. Either may serve as a backend implementation.
- **Adopting TVM or MLIR wholesale.** Rejected as the contract. Both are large dependencies that keep changing, and TVM is oriented towards tensor workloads in machine learning. Either may be used inside a backend.
- **Plain C++ only.** Rejected: the only path to a new device would be rewriting.

## 4. Consequences

- **The engineering effort is substantial.** The stack must start small (an IR, a CPU interpreter and one kernel) and grow with the vertical slices.
- **Every backend can be tested** against the reference interpreter.
- **A second backend in Phase 2** tests whether the abstraction really holds.
- **Units and frames could become IR types.** The IR is a natural place to type them ([idea](../ideas/2026-10-21-units-and-frames-as-types.md)).
